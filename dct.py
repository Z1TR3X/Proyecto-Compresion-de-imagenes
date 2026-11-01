"""
dct.py  -  Compresión de imágenes con la DCT global (base ortonormal).

Qué hace:
  1. Conserva 1, 5, 10 y 25 % de los coeficientes y muestra las reconstrucciones.
  2. Compara dos reglas: soporte PREFIJADO (frecuencias bajas) y selección por MAGNITUD.
  3. Calcula error relativo, MSE, PSNR y energía retenida (antes de cuantizar).
  4. Cuantiza los coeficientes y estima los bits (valores + posiciones + cabecera).

Ejecutar:  python dct.py
Lee las imágenes de la carpeta  imagenes/  y guarda todo en  resultados/
"""
import os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")                       # guarda las figuras sin abrir ventanas
import matplotlib.pyplot as plt
from PIL import Image
from scipy.fft import dctn, idctn

TAM = 256                       # tamaño de trabajo: 256 x 256
PORCENTAJES = [1, 5, 10, 25]    # % de coeficientes que se conservan
DELTA = 0.05                    # paso de cuantización
IMAGEN = ""                     # "" = usa todas las imágenes de imagenes/; o pon un nombre, ej. "radiografia.png"


# ---------- 1. Cargar imagen: gris, cuadrada, 256x256, valores en [0, 1] ----------
def cargar_imagen(ruta):
    img = Image.open(ruta).convert("L")
    w, h = img.size
    lado = min(w, h)
    img = img.crop(((w - lado) // 2, (h - lado) // 2, (w + lado) // 2, (h + lado) // 2))
    img = img.resize((TAM, TAM), Image.LANCZOS)
    return np.asarray(img, dtype=float) / 255.0


# ---------- 2. DCT ortonormal ----------
# Matriz D: solo se usa para VERIFICAR con una matriz pequeña (sección 6, paso 2 del documento)
def matriz_dct(N):
    """D[k, l] = a_k * cos(pi * k * (2l + 1) / (2N)),  a_0 = sqrt(1/N),  a_k = sqrt(2/N)"""
    k = np.arange(N).reshape(-1, 1)
    l = np.arange(N).reshape(1, -1)
    D = np.sqrt(2 / N) * np.cos(np.pi * k * (2 * l + 1) / (2 * N))
    D[0, :] = np.sqrt(1 / N)
    return D


# Transformada rápida de SciPy (da lo mismo que C = D A D^T, pero más rápido)
def transformar(A):  return dctn(A, norm="ortho")     # C = D A D^T
def reconstruir(C):  return idctn(C, norm="ortho")    # A = D^T C D


# ---------- 3. Reglas de selección (devuelven una máscara: True = se conserva) ----------
def seleccion_magnitud(C, k):
    """Los k coeficientes de mayor |c|."""
    orden = np.argsort(np.abs(C).ravel())[-k:]     # posiciones de los k más grandes
    M = np.zeros(C.size, dtype=bool)
    M[orden] = True
    return M.reshape(C.shape)


def seleccion_prefijada(C, k):
    """Las k frecuencias más bajas (esquina superior izquierda, ordenadas por i + j).
       Las posiciones no dependen de la imagen."""
    i, j = np.indices(C.shape)
    orden = np.argsort((i + j).ravel(), kind="stable")[:k]
    M = np.zeros(C.size, dtype=bool)
    M[orden] = True
    return M.reshape(C.shape)


# ---------- 4. Métricas y bits ----------
def metricas(A, A_rec):
    mse = np.mean((A - A_rec) ** 2)
    error_rel = np.linalg.norm(A - A_rec) / np.linalg.norm(A)
    psnr = np.inf if mse == 0 else 10 * np.log10(1 / mse)
    return error_rel, mse, psnr


def bits_estimados(q, M, regla):
    """cabecera (64) + valores (k enteros) + posiciones (solo si la regla es por magnitud)."""
    k = int(M.sum())
    bits_por_valor = int(np.ceil(np.log2(np.abs(q[M]).max() + 1))) + 1   # +1 por el signo
    bits_posiciones = k * int(np.log2(q.size)) if regla == "magnitud" else 0   # bits de una posición (16 en 256x256)
    return 64 + k * bits_por_valor + bits_posiciones


# ---------- 5. Verificaciones ----------
print("Verificaciones DCT")
D8 = matriz_dct(8)                          # matriz pequeña de 8 x 8
np.random.seed(0)                           # semilla fija: X sale igual en cada ejecución
X = np.random.rand(8, 8)
print("  D D^T = I  (8x8)            ->", np.allclose(D8 @ D8.T, np.eye(8)))
print("  dctn(X) = D X D^T  (8x8)    ->", np.allclose(transformar(X), D8 @ X @ D8.T))

# ---------- 6. Experimento ----------
os.makedirs("resultados", exist_ok=True)
EXTENSIONES = (".png", ".jpg", ".jpeg", ".webp", ".bmp", ".tif", ".tiff")
rutas = sorted(r for r in glob.glob("imagenes/*") if r.lower().endswith(EXTENSIONES))
if IMAGEN:
    rutas = [os.path.join("imagenes", IMAGEN)]     # solo la imagen elegida

for ruta in rutas:
    nombre = os.path.splitext(os.path.basename(ruta))[0]
    A = cargar_imagen(ruta)
    C = transformar(A)
    print(f"\nImagen: {nombre}")
    print("  reconstrucción exacta     ->", np.allclose(reconstruir(C), A))
    print("  Parseval ||A||² = ||C||²  ->", np.isclose(np.sum(A**2), np.sum(C**2)))
    print(f"  Tamaño sin comprimir: {A.size * 8 / 8 / 1024:.1f} KB  (8 bits por píxel)")
    print(f"  {'regla':<10}{'%':>4}{'err_rel':>10}{'PSNR':>9}{'energía':>10}{'err=desc':>10}{'PSNR_q':>9}{'bpp':>7}{'KB':>8}")

    fig, ejes = plt.subplots(2, 5, figsize=(16, 7))
    for f, (regla, seleccionar) in enumerate([("prefijada", seleccion_prefijada), ("magnitud", seleccion_magnitud)]):
        ejes[f, 0].imshow(A, cmap="gray", vmin=0, vmax=1)
        ejes[f, 0].set_title("Original")
        ejes[f, 0].set_ylabel(regla, fontsize=12)

        for c, p in enumerate(PORCENTAJES):
            k = round(A.size * p / 100)
            M = seleccionar(C, k)
            A_rec = reconstruir(C * M)

            # Métricas antes de cuantizar
            error_rel, mse, psnr = metricas(A, A_rec)
            energia = np.sum(C[M] ** 2) / np.sum(C ** 2)
            err_igual_desc = np.isclose(np.sum((A - A_rec) ** 2), np.sum(C[~M] ** 2))

            # Cuantización: q = round(c / Delta),  c_hat = Delta * q
            q = np.round(C / DELTA) * M
            _, _, psnr_q = metricas(A, reconstruir(DELTA * q))
            bits = bits_estimados(q, M, regla)
            bpp = bits / A.size                               # bits por píxel (sin comprimir = 8)
            kb = bits / 8 / 1024                              # tamaño estimado comprimido en KB

            print(f"  {regla:<10}{p:>4}{error_rel:>10.4f}{psnr:>9.2f}{energia:>10.5f}{str(err_igual_desc):>10}{psnr_q:>9.2f}{bpp:>7.3f}{kb:>8.1f}")

            ejes[f, c + 1].imshow(np.clip(A_rec, 0, 1), cmap="gray", vmin=0, vmax=1)   # misma escala
            ejes[f, c + 1].set_title(f"{p}%  PSNR={psnr:.1f} dB")

    for ax in ejes.ravel():
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f"DCT - {nombre}")
    fig.tight_layout()
    fig.savefig(f"resultados/dct_{nombre}.png", dpi=100)
    plt.close(fig)

if rutas:
    print("\nListo: figuras en resultados/")
else:
    print("No hay imágenes en la carpeta imagenes/")
