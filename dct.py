MOSTRAR_SIN_DELTA = False 
"""
dct.py  -  Compresión de imágenes con la DCT global (base ortonormal).
"""
import os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")                       # guarda las figuras sin abrir ventanas
import matplotlib.pyplot as plt
from PIL import Image
from scipy.fft import dctn, idctn

PORCENTAJES = [1, 5, 10, 25]    
DELTA = 0.05                    
RUTA = "imagenes/texturayperiodicidad.jpeg"   
GUARDAR_PREFIJADA = True    
GUARDAR_MAGNITUD = True         
GUARDAR_CUADRO = True           
MOSTRAR_SIN_DELTA = False       


# ----- 1. Cargar imagen: gris, valores en [0, 1] ------
def cargar_imagen(ruta):
    img = Image.open(ruta).convert("L")
    return np.asarray(img, dtype=float) / 255.0


# --- 2. DCT ortonormal ------
# Matriz D:para VERIFICAR 
def matriz_dct(N):
    """D[k, l] = a_k * cos(pi * k * (2l + 1) / (2N)),  a_0 = sqrt(1/N),  a_k = sqrt(2/N)"""
    k = np.arange(N).reshape(-1, 1)
    l = np.arange(N).reshape(1, -1)
    D = np.sqrt(2 / N) * np.cos(np.pi * k * (2 * l + 1) / (2 * N))
    D[0, :] = np.sqrt(1 / N)
    return D


# Transformada rápida de SciPy 
def transformar(A):  return dctn(A, norm="ortho")     # C = D A D^T
def reconstruir(C):  return idctn(C, norm="ortho")    # A = D^T C D


# ------- 3. Reglas de selección  ---
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


# -4. Error de recuperación--
def error_recuperacion(A, A_rec):
    """||A - A_rec||_F / ||A||_F"""
    return np.linalg.norm(A - A_rec) / np.linalg.norm(A)


def bits_estimados(q, M, regla):
    """cabecera (64) + valores (k enteros) + posiciones (solo si la regla es por magnitud)."""
    k = int(M.sum())
    bits_por_valor = int(np.ceil(np.log2(np.abs(q[M]).max() + 1))) + 1   # +1 por el signo
    bits_posiciones = k * int(np.ceil(np.log2(q.size))) if regla == "magnitud" else 0   # bits de una posición (redondeado hacia arriba)
    return 64 + k * bits_por_valor + bits_posiciones


# ------ 5. Verificaciones -
print("Verificaciones DCT")
D8 = matriz_dct(8)                          
np.random.seed(0)                           
X = np.random.rand(8, 8)
print("  D D^T = I  (8x8)            ->", np.allclose(D8 @ D8.T, np.eye(8)))
print("  dctn(X) = D X D^T  (8x8)    ->", np.allclose(transformar(X), D8 @ X @ D8.T))

# ------ 6. Experimento ----
os.makedirs("resultados", exist_ok=True)
rutas = [RUTA]                  # solo la imagen elegida

for ruta in rutas:
    nombre = os.path.splitext(os.path.basename(ruta))[0]
    A = cargar_imagen(ruta)
    C = transformar(A)
    print(f"\nImagen: {nombre}")
    print("  Parseval ||A||² = ||C||²  ->", np.isclose(np.sum(A**2), np.sum(C**2)))
    print(f"  Tamaño real: {A.shape[0]} x {A.shape[1]} píxeles = {A.size / 1024:.1f} KB sin comprimir  (8 bits por píxel)")
    extra = f"{'error sin Δ':>13}{'MSE sin Δ':>12}{'energía retenida':>18}" if MOSTRAR_SIN_DELTA else ""
    print(f"  {'regla':<10}{'% coef':>7}{extra}{'MSE':>12}{'error de recuperación':>24}{'tasa de compresión':>21}{'tamaño comprimido':>20}")

    fig, ejes = plt.subplots(2, 5, figsize=(16, 6.4 * A.shape[0] / A.shape[1] + 1.2))   # alto según la proporción de la imagen
    for f, (regla, seleccionar) in enumerate([("prefijada", seleccion_prefijada), ("magnitud", seleccion_magnitud)]):
        ejes[f, 0].imshow(A, cmap="gray", vmin=0, vmax=1)
        ejes[f, 0].set_title("Original")
        ejes[f, 0].set_ylabel(regla, fontsize=12)

        for c, p in enumerate(PORCENTAJES):
            k = round(A.size * p / 100)
            M = seleccionar(C, k)

            # Cuantización: q = round(c / Delta),  c_hat = Delta * q
            q = np.round(C / DELTA) * M
            A_rec = reconstruir(DELTA * q)                    # imagen recuperada a partir de lo que se guarda

            A_sin = reconstruir(C * M)                               # sin Δ: solo descarte de coeficientes
            error_sin = error_recuperacion(A, A_sin)
            mse_sin = np.mean((A - A_sin) ** 2)                      # MSE sin Δ
            energia = np.sum(C[M] ** 2) / np.sum(C ** 2)             # energía retenida
            error_rel = error_recuperacion(A, A_rec)
            mse = np.mean((A - A_rec) ** 2)                          # MSE: promedio por píxel de (A - A_rec)^2
            bits = bits_estimados(q, M, regla)
            tasa = A.size * 8 / bits                          # tasa de compresión: bits originales / bits comprimidos
            kb = bits / 8 / 1024                              # tamaño estimado comprimido en KB

            extra = f"{error_sin:>13.2%}{mse_sin:>12.6f}{energia:>18.2%}" if MOSTRAR_SIN_DELTA else ""
            print(f"  {regla:<10}{p:>7}{extra}{mse:>12.6f}{error_rel:>24.2%}{tasa:>19.1f}:1{kb:>17.1f} KB")

            # guardar cada imagen comprimida (tamaño real) en resultados/
            if (regla == "prefijada" and GUARDAR_PREFIJADA) or (regla == "magnitud" and GUARDAR_MAGNITUD):
                Image.fromarray(np.uint8(np.round(np.clip(A_rec, 0, 1) * 255))).save(f"resultados/{nombre}_{regla}_{p}.png")

            ejes[f, c + 1].imshow(np.clip(A_rec, 0, 1), cmap="gray", vmin=0, vmax=1)   # misma escala
            ejes[f, c + 1].set_title(f"{p}% de coeficientes")

    for ax in ejes.ravel():
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f"DCT - {nombre}")
    fig.tight_layout()
    if GUARDAR_CUADRO:
        fig.savefig(f"resultados/dct_{nombre}.png", dpi=100)
    plt.close(fig)

if rutas:
    print("\nListo: figuras en resultados/")
else:
    print("No hay imágenes en la carpeta imagenes/")
