import os, glob, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

TAM_BLOQUE = 8
PORCENTAJES = [1, 5, 10, 15, 25]     # 1, 5, 10, 25 % del grupo + 15 % que usa dct_escala_real.py
CARPETA = "resultados_escala_real"
NOMBRES = {"Gradiente": "Resonancia", "bordesmarcados": "Radiografía", "texturayperiodicidad": "Tomografía"}

BITS_VALOR = 8                                  # bits por coeficiente conservado (igual que dct_escala_real: bpp = 8 / CR)
BITS_POSICION = int(np.log2(TAM_BLOQUE ** 2))   # 6 bits para ubicar un coeficiente dentro del bloque (solo en "magnitud")


# ---------------Matriz ortonormal de Haar H_N (recursiva) ----------
def etapa_haar(n):
    """Una etapa: primeras n/2 filas = promedios (x1 + x2)/sqrt(2); últimas n/2 = diferencias (x1 - x2)/sqrt(2)."""
    E = np.zeros((n, n))
    for i in range(n // 2):
        E[i, 2*i] = E[i, 2*i + 1] = 1 / np.sqrt(2)
        E[n//2 + i, 2*i], E[n//2 + i, 2*i + 1] = 1 / np.sqrt(2), -1 / np.sqrt(2)
    return E


def construir_matriz_haar(N):
    """H_N = diag(H_{N/2}, I) @ etapa(N), con H_1 = [1]. Para N = 8 es la matriz H_8 del informe."""
    if N == 1:
        return np.ones((1, 1))
    H = np.eye(N)
    H[:N//2, :N//2] = construir_matriz_haar(N // 2)
    return H @ etapa_haar(N)


# ---------- Reglas de selección ----------
def mascara_prefijada(tam, k):
    """Primeros k coeficientes en orden de cuadrados crecientes (escala 1, 2, 4, 8) en la esquina superior izquierda."""
    i, j = np.indices((tam, tam))
    orden = np.lexsort(((i + j).ravel(), np.maximum(i, j).ravel()))[:k]
    m = np.zeros(tam * tam, dtype=bool)
    m[orden] = True
    return m.reshape(tam, tam)


def mascara_magnitud(bloques_C, k):
    """Para cada bloque, True en sus k coeficientes de mayor |c|."""
    plano = np.abs(bloques_C).reshape(*bloques_C.shape[:2], -1)
    idx = np.argsort(plano, axis=-1)[..., -k:]
    m = np.zeros_like(plano, dtype=bool)
    np.put_along_axis(m, idx, True, axis=-1)
    return m.reshape(bloques_C.shape)


# --------- Imagen <-> bloques ----------
def a_bloques(U, b):
    M, N = U.shape
    return U.reshape(M // b, b, N // b, b).swapaxes(1, 2)          # (M/b, N/b, b, b)


def de_bloques(B):
    nb_i, nb_j, b, _ = B.shape
    return B.swapaxes(1, 2).reshape(nb_i * b, nb_j * b)


# --------Verificaciones ----------
H = construir_matriz_haar(TAM_BLOQUE)
HT = H.T
print("Verificaciones de la base de Haar")
print(f"  H_{TAM_BLOQUE} H_{TAM_BLOQUE}^T = I   ->", np.allclose(H @ HT, np.eye(TAM_BLOQUE)))


# ----------Experimento ----------
os.makedirs(CARPETA, exist_ok=True)
EXTENSIONES = (".png", ".jpg", ".jpeg", ".webp", ".bmp", ".tif", ".tiff")
rutas = sorted(r for r in glob.glob("imagenes/*") if r.lower().endswith(EXTENSIONES))
filas_tabla = []

for ruta in rutas:
    archivo = os.path.splitext(os.path.basename(ruta))[0]
    nombre = NOMBRES.get(archivo, archivo)
    peso_archivo_kb = os.path.getsize(ruta) / 1024
    U = np.array(Image.open(ruta).convert("L"), dtype=np.float64)       # 0..255, tamaño real
    M, N = U.shape
    pad_m = (TAM_BLOQUE - M % TAM_BLOQUE) % TAM_BLOQUE
    pad_n = (TAM_BLOQUE - N % TAM_BLOQUE) % TAM_BLOQUE
    U_pad = np.pad(U, ((0, pad_m), (0, pad_n)), mode="edge")
    M_pad, N_pad = U_pad.shape
    n_bloques = (M_pad // TAM_BLOQUE) * (N_pad // TAM_BLOQUE)

    B = a_bloques(U_pad, TAM_BLOQUE)
    C = H @ B @ HT                                                       # C = H A H^T en cada bloque
    parseval = np.isclose(np.sum(B ** 2), np.sum(C ** 2))
    exacta = np.allclose(de_bloques(HT @ C @ H)[:M, :N], U)

    bytes_originales = M * N
    print(f"\nImagen: {nombre}  ({archivo})")
    print(f"  Tamaño real            : {N} x {M} píxeles (ancho x alto);  M = {M} filas, N = {N} columnas")
    print(f"  Con padding edge       : {N_pad} x {M_pad} píxeles, {n_bloques} bloques de 8x8 (se recorta al tamaño real al reconstruir)")
    print(f"  Peso original (M*N)    : {bytes_originales / 1024:.1f} KB  ({bytes_originales:,} bytes, 8 bits por píxel)")
    print(f"  Archivo en carpeta     : {peso_archivo_kb:.1f} KB")
    print(f"  Parseval ||A||_F² = ||C||_F² (por bloques)  ->", parseval)
    print(f"  Reconstrucción exacta al 100 % (recortada)  ->", exacta)

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
            U_rec = np.clip(de_bloques(HT @ Mc @ H)[:M, :N], 0, 255)      # A_rec = H^T M H, recorte al tamaño original
            err = np.linalg.norm(U - U_rec) / np.linalg.norm(U)
            mse = np.mean((U - U_rec) ** 2)
            coef_totales = k * n_bloques
            bits = coef_totales * (BITS_VALOR + (BITS_POSICION if regla == "magnitud" else 0))
            kb, bpp = bits / 8 / 1024, bits / (M * N)
            filas_tabla.append((nombre, f"{N}x{M}", regla, p, k, err, mse, kb, bpp, bytes_originales / 1024))
            ejes[f, c + 1].imshow(U_rec, cmap="gray", vmin=0, vmax=255)
            ejes[f, c + 1].set_title(f"{p}% | {err:.4f} | {mse:.2f}", fontsize=10)
    for ax in ejes.ravel():
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f"Haar 8x8 por bloques - {nombre} ({N}x{M})   (subtítulos: % Coef | Error Rel | MSE)", fontsize=11)
    fig.tight_layout(h_pad=1.5)
    fig.savefig(f"{CARPETA}/haar_escala_real_{archivo}.png", dpi=90)
    plt.close(fig)

# ----------Tabla en consola, CSV y Markdown ----------
cab = ["Imagen", "Regla", "% Coef", "Error Relativo", "MSE", "Tamaño estimado (KB)"]
print(f"\nTamaño estimado = k·bloques·({BITS_VALOR} bits por valor" + f" + {BITS_POSICION} bits de posición solo en 'magnitud') / 8 / 1024.  k = coeficientes por bloque (de 64).")
print(f"{cab[0]:<12}| {cab[1]:<10}| {cab[2]:>11} | {cab[3]:>14} | {cab[4]:>9} | {cab[5]:>20} | {'Peso original (KB)':>18}")
print("-" * 12 + "+" + "-" * 11 + "+" + "-" * 13 + "+" + "-" * 16 + "+" + "-" * 11 + "+" + "-" * 22 + "+" + "-" * 19)
for nombre, dim, regla, p, k, err, mse, kb, bpp, orig in filas_tabla:
    print(f"{nombre:<12}| {regla:<10}| {f'{p}% (k={k})':>11} | {err:>14.4f} | {mse:>9.2f} | {kb:>20.1f} | {orig:>18.1f}")

with open(f"{CARPETA}/tabla_haar_escala_real.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["Imagen", "Dimensiones", "Regla", "% Coef (nominal)", "k por bloque", "% real", "Error Relativo", "MSE", "Tamaño estimado (KB)", "bpp", "Peso original (KB)"])
    for nombre, dim, regla, p, k, err, mse, kb, bpp, orig in filas_tabla:
        w.writerow([nombre, dim, regla, p, k, f"{100 * k / TAM_BLOQUE**2:.2f}", f"{err:.4f}", f"{mse:.2f}", f"{kb:.1f}", f"{bpp:.3f}", f"{orig:.1f}"])
with open(f"{CARPETA}/tabla_haar_escala_real.md", "w", encoding="utf-8") as fh:
    fh.write("| " + " | ".join(cab) + " | bpp | Peso original (KB) |\n|" + "---|" * (len(cab) + 2) + "\n")
    for nombre, dim, regla, p, k, err, mse, kb, bpp, orig in filas_tabla:
        fh.write(f"| {nombre} ({dim}) | {regla} | {p}% (k={k}) | {err:.4f} | {mse:.2f} | {kb:.1f} | {bpp:.3f} | {orig:.1f} |\n")

print(f"\nListo: figuras y tablas en {CARPETA}/" if rutas else "No hay imágenes en la carpeta imagenes/")
