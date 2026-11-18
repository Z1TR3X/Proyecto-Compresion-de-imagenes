import os
import glob
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

TAM_BLOQUE = 8
PORCENTAJES = [1, 5, 10, 25]
CARPETA = "resultados_escala_real"
NOMBRES = {
    "Gradiente": "Resonancia",
    "bordesmarcados": "Radiografía",
    "texturayperiodicidad": "Tomografía",
}

BITS_VALOR = 8                                  
BITS_POSICION = int(np.log2(TAM_BLOQUE ** 2))   


# ---------- Matriz ortonormal de Haar H_N (explícita, recursiva) ----------
def etapa_haar(n: int) -> np.ndarray:
    """Promedios (x1 + x2)/sqrt(2) en la mitad superior y diferencias (x1 - x2)/sqrt(2) en la inferior."""
    E = np.zeros((n, n), dtype=np.float64)
    inv_sqrt2 = 1.0 / np.sqrt(2.0)
    for i in range(n // 2):
        E[i, 2 * i] = inv_sqrt2
        E[i, 2 * i + 1] = inv_sqrt2
        E[n // 2 + i, 2 * i] = inv_sqrt2
        E[n // 2 + i, 2 * i + 1] = -inv_sqrt2
    return E


def construir_matriz_haar(n: int) -> np.ndarray:
    """Construye H_N = diag(H_{N/2}, I_{N/2}) @ etapa_haar(N), con H_1 = [1]."""
    if n == 1:
        return np.ones((1, 1), dtype=np.float64)
    H = np.eye(n, dtype=np.float64)
    H[: n // 2, : n // 2] = construir_matriz_haar(n // 2)
    return H @ etapa_haar(n)


# ---------- Reglas de selección por bloque ----------
def mascara_prefijada(tam: int, k: int) -> np.ndarray:
    """Primeros k coeficientes en esquina superior izquierda (escalas 1, 2, 4, 8)."""
    i, j = np.indices((tam, tam))
    orden = np.lexsort(((i + j).ravel(), np.maximum(i, j).ravel()))[:k]
    m = np.zeros(tam * tam, dtype=bool)
    m[orden] = True
    return m.reshape(tam, tam)


def mascara_magnitud(bloques_C: np.ndarray, k: int) -> np.ndarray:
    """Para cada bloque 8x8, selecciona los k coeficientes de mayor valor absoluto |c|."""
    plano = np.abs(bloques_C).reshape(*bloques_C.shape[:2], -1)
    idx = np.argsort(plano, axis=-1)[..., -k:]
    m = np.zeros_like(plano, dtype=bool)
    np.put_along_axis(m, idx, True, axis=-1)
    return m.reshape(bloques_C.shape)


# ---------- Conversión imagen <-> bloques ----------
def a_bloques(U: np.ndarray, b: int) -> np.ndarray:
    M, N = U.shape
    return U.reshape(M // b, b, N // b, b).swapaxes(1, 2)


def de_bloques(B: np.ndarray) -> np.ndarray:
    nb_i, nb_j, b, _ = B.shape
    return B.swapaxes(1, 2).reshape(nb_i * b, nb_j * b)


# ---------- Verificación algebraica de la base ----------
H = construir_matriz_haar(TAM_BLOQUE)
HT = H.T
assert np.allclose(H @ HT, np.eye(TAM_BLOQUE)), "H_8 no es ortonormal."
print(f"Verificación H_{TAM_BLOQUE} @ H_{TAM_BLOQUE}^T = I_{TAM_BLOQUE} -> True")


# ---------- Procesamiento en dimensión real ----------
os.makedirs(CARPETA, exist_ok=True)
EXTENSIONES = (".png", ".jpg", ".jpeg", ".webp", ".bmp", ".tif", ".tiff")
rutas = sorted(r for r in glob.glob("imagenes/*") if r.lower().endswith(EXTENSIONES))
filas_tabla = []

for ruta in rutas:
    archivo = os.path.splitext(os.path.basename(ruta))[0]
    nombre = NOMBRES.get(archivo, archivo)
    U = np.array(Image.open(ruta).convert("L"), dtype=np.float64)
    M, N = U.shape

    pad_m = (TAM_BLOQUE - M % TAM_BLOQUE) % TAM_BLOQUE
    pad_n = (TAM_BLOQUE - N % TAM_BLOQUE) % TAM_BLOQUE
    U_pad = np.pad(U, ((0, pad_m), (0, pad_n)), mode="edge")
    M_pad, N_pad = U_pad.shape
    n_bloques = (M_pad // TAM_BLOQUE) * (N_pad // TAM_BLOQUE)

    B = a_bloques(U_pad, TAM_BLOQUE)
    C = H @ B @ HT
    parseval = np.isclose(np.sum(B ** 2), np.sum(C ** 2))
    exacta = np.allclose(de_bloques(HT @ C @ H)[:M, :N], U)

    bytes_originales = M * N
    print(f"\nImagen: {nombre} ({N}x{M} px | Peso original: {bytes_originales / 1024:.1f} KB)")
    print(f"  Parseval ||A||_F^2 = ||C||_F^2 (por bloques) -> {parseval}")
    print(f"  Reconstrucción exacta al 100% (recortada)    -> {exacta}")

    ncol = len(PORCENTAJES) + 1
    fig, ejes = plt.subplots(2, ncol, figsize=(2.6 * ncol, 2 * 2.6 * M / N + 0.9))

    for f, regla in enumerate(["prefijada", "magnitud"]):
        ejes[f, 0].imshow(U, cmap="gray", vmin=0, vmax=255)
        ejes[f, 0].set_title("Original", fontsize=10)
        ejes[f, 0].set_ylabel(regla, fontsize=11)

        for c, p in enumerate(PORCENTAJES):
            k = max(1, int(TAM_BLOQUE ** 2 * p / 100))
            if regla == "prefijada":
                masc = np.broadcast_to(mascara_prefijada(TAM_BLOQUE, k), C.shape)
            else:
                masc = mascara_magnitud(C, k)

            Mc = C * masc
            U_rec = np.clip(de_bloques(HT @ Mc @ H)[:M, :N], 0, 255)

            err = np.linalg.norm(U - U_rec, ord="fro") / np.linalg.norm(U, ord="fro")
            mse = np.mean((U - U_rec) ** 2)

            coef_totales = k * n_bloques
            bits = coef_totales * (BITS_VALOR + (BITS_POSICION if regla == "magnitud" else 0))
            kb = bits / 8 / 1024
            bpp = bits / (M * N)

            filas_tabla.append((nombre, f"{N}x{M}", regla, p, k, err, mse, kb, bpp, bytes_originales / 1024))

            ejes[f, c + 1].imshow(U_rec, cmap="gray", vmin=0, vmax=255)
            ejes[f, c + 1].set_title(f"{p}% | {err:.4f} | {mse:.2f}", fontsize=10)

    for ax in ejes.ravel():
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle(f"Haar 8x8 por bloques - {nombre} ({N}x{M})  (% Coef | Error Rel | MSE)", fontsize=11)
    fig.tight_layout(h_pad=1.5)
    fig.savefig(f"{CARPETA}/haar_escala_real_{archivo}.png", dpi=110)
    plt.close(fig)

# ---------- Tabla limpia en consola (solo 6 columnas) y CSV ----------
print("\n" + "=" * 86)
print(f"{'Imagen':<12}| {'Regla':<10}| {'% Coeficientes':>15} | {'Error Relativo':>14} | {'MSE':>9} | {'Tamaño estimado (KB)':>20}")
print("-" * 86)
for nombre, dim, regla, p, k, err, mse, kb, bpp, orig in filas_tabla:
    print(f"{nombre:<12}| {regla:<10}| {f'{p}% (k={k})':>15} | {err:>14.4f} | {mse:>9.2f} | {kb:>20.1f}")
print("=" * 86)

with open(f"{CARPETA}/tabla_haar_escala_real.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["Imagen", "Dimensiones", "Regla", "% Coeficientes", "k por bloque", "Error Relativo", "MSE", "Tamaño estimado (KB)", "bpp", "Peso original (KB)"])
    for nombre, dim, regla, p, k, err, mse, kb, bpp, orig in filas_tabla:
        w.writerow([nombre, dim, regla, f"{p}%", k, f"{err:.4f}", f"{mse:.2f}", f"{kb:.1f}", f"{bpp:.3f}", f"{orig:.1f}"])