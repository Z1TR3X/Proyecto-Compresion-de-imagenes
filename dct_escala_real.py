import os
import numpy as np
from PIL import Image

def construir_matriz_dct(N: int) -> np.ndarray:
    """Construye la matriz ortogonal C de tamaño N x N de la DCT-II."""
    C = np.zeros((N, N), dtype=np.float64)
    for k in range(N):
        alpha = np.sqrt(1.0 / N) if k == 0 else np.sqrt(2.0 / N)
        for n in range(N):
            C[k, n] = alpha * np.cos(np.pi * (2 * n + 1) * k / (2 * N))
    return C

def comprimir_imagen_dct(
    ruta_imagen: str, tam_bloque: int = 8, porcentaje_conservar: float = 0.15
):
    """Comprime y reconstruye una imagen en escala de grises manteniendo su tamaño M x N."""
    # 1. Cargar imagen en escala de grises (8 bits)
    peso_archivo_entrada_kb = os.path.getsize(ruta_imagen) / 1024
    img = Image.open(ruta_imagen).convert("L")
    U_orig = np.array(img, dtype=np.float64)
    M, N = U_orig.shape

    # Padding si M o N no son múltiplos del tamaño de bloque (para conservar tamaño exacto)
    pad_m = (tam_bloque - M % tam_bloque) % tam_bloque
    pad_n = (tam_bloque - N % tam_bloque) % tam_bloque
    U_pad = np.pad(U_orig, ((0, pad_m), (0, pad_n)), mode="edge")
    M_pad, N_pad = U_pad.shape

    # 2. Construir matriz ortogonal real C (C^-1 = C^T)
    C = construir_matriz_dct(tam_bloque)
    C_T = C.T

    # Número de coeficientes a retener por bloque
    total_coef_bloque = tam_bloque * tam_bloque
    k_retener = max(1, int(total_coef_bloque * porcentaje_conservar))

    U_V_comprimida = np.zeros_like(U_pad)
    U_rec_pad = np.zeros_like(U_pad)
    coeficientes_no_cero = 0

    # 3. Procesamiento por bloques N_b x N_b
    for i in range(0, M_pad, tam_bloque):
        for j in range(0, N_pad, tam_bloque):
            bloque_U = U_pad[i : i + tam_bloque, j : j + tam_bloque]

            # Transformada Directa 2D separable: V = C * U * C^T
            bloque_V = C @ bloque_U @ C_T

            # Umbralización (descarta altas frecuencias de baja energía)
            V_flat = np.abs(bloque_V.flatten())
            umbral = np.partition(V_flat, -k_retener)[-k_retener]
            bloque_V_comp = np.where(np.abs(bloque_V) >= umbral, bloque_V, 0.0)

            U_V_comprimida[i : i + tam_bloque, j : j + tam_bloque] = (
                bloque_V_comp
            )
            coeficientes_no_cero += np.count_nonzero(bloque_V_comp)

            # Transformada Inversa 2D: U_hat = C^T * V_hat * C
            bloque_U_rec = C_T @ bloque_V_comp @ C
            U_rec_pad[i : i + tam_bloque, j : j + tam_bloque] = bloque_U_rec

    # 4. Guardar el archivo comprimido real (coeficientes dispersos)
    V_enteros = np.round(U_V_comprimida).astype(np.int16)
    np.savez_compressed("archivo_comprimido_real.npz", dct=V_enteros)
    peso_comprimido_npz_kb = (
        os.path.getsize("archivo_comprimido_real.npz") / 1024
    )

    # 5. Recortar padding y guardar imagen reconstruida para visualización
    U_rec = np.clip(U_rec_pad[:M, :N], 0, 255)
    Image.fromarray(U_rec.astype(np.uint8)).save(
        "imagen_reconstruida_dct.jpg", quality=75
    )
    peso_reconstruido_jpg_kb = (
        os.path.getsize("imagen_reconstruida_dct.jpg") / 1024
    )

    # 6. Cálculo de Fidelidad (MSE y PSNR) y Tasas de Compresión
    mse = np.mean((U_orig - U_rec) ** 2)
    psnr = 10.0 * np.log10((255.0**2) / mse) if mse > 0 else float("inf")

    total_coef = M_pad * N_pad
    tasa_compresion_cr = total_coef / coeficientes_no_cero
    porcentaje_reduccion = (1.0 - (coeficientes_no_cero / total_coef)) * 100.0
    bpp = 8.0 / tasa_compresion_cr

    bytes_originales = M * N
    bytes_comprimidos = bytes_originales / tasa_compresion_cr

    print("--- TAMAÑOS EN DISCO Y MEMORIA ---")
    print(
        f"Archivo original en carpeta (.jpg)   : {peso_archivo_entrada_kb:.2f} KB"
    )
    print(
        f"Imagen reconstruida en carpeta (.jpg): {peso_reconstruido_jpg_kb:.2f} KB"
    )
    print(
        f"Archivo comprimido real (.npz)       : {peso_comprimido_npz_kb:.2f} KB"
    )
    print(
        f"Matriz cruda sin comprimir (RAM)     : {bytes_originales / 1024:.2f} KB ({bytes_originales:,} bytes)"
    )
    print(
        f"Flujo teórico comprimido DCT         : {bytes_comprimidos / 1024:.2f} KB ({bytes_comprimidos:,.0f} bytes)\n"
    )

    print("--- MÉTRICAS PARA EL INFORME ---")
    print(f"Dimensiones originales  : {M} x {N} píxeles")
    print(f"Dimensiones finales     : {U_rec.shape[0]} x {U_rec.shape[1]} píxeles")
    print(f"Coeficientes eliminados : {porcentaje_reduccion:.2f}%")
    print(f"Tasa de Compresión (CR) : {tasa_compresion_cr:.2f} : 1")
    print(f"Bits por píxel (bpp)    : {bpp:.2f} bpp (de 8.00 bpp originales)")
    print(f"MSE (sigma_ls^2)        : {mse:.4f}")
    print(f"PSNR (SNR')             : {psnr:.2f} dB")


comprimir_imagen_dct(
    "bordesmarcados.jpg", tam_bloque=8, porcentaje_conservar=0.15
)