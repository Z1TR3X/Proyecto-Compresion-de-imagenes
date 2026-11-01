# Compresión de radiografías con bases ortonormales: DCT vs Haar

**Curso:** Álgebra Lineal Aplicada 2026-2 — Universidad del Pacífico
**Tema:** Compresión de imágenes mediante bases ortogonales


## Problema y objetivo
comparar qué base ortonormal (DCT global o Haar multinivel) conserva mejor la información de distintos tipos de radiografía (suave, con bordes y con textura) cuando se guarda solo una fracción de los coeficientes, medido con error relativo, PSNR, energía retenida y bits estimados.


## Organización de archivos
```
dct.py                -> experimento con la DCT global
haar.py               -> experimento con Haar multinivel
requirements.txt      -> dependencias
imagenes/             -> imágenes de entrada (PNG/JPG)
resultados/           -> figuras que generan dct.py y haar.py
```

## Cómo ejecutarlo
Requisitos: Python 3.10 o superior.

```bash
pip install -r requirements.txt
python dct.py                   # experimento con DCT
python haar.py                  # experimento con Haar
```

Coloquen sus imágenes (PNG preferido) en `imagenes/`. El script las convierte a escala de grises, recorta un cuadrado central, las lleva a 256×256 (potencia de 2, necesario para Haar) y normaliza al rango [0, 1].

## Qué hacen `dct.py` y `haar.py`
Los dos archivos siguen los mismos pasos, cada uno con su base: la DCT [1], [2] y Haar ortonormal multinivel [3]. El procedimiento sigue el material del curso [4].

1. **Verificaciones:** ortogonalidad de la matriz, reconstrucción exacta y Parseval en cada imagen. `dct.py` calcula la DCT con `scipy.fft.dctn` [5] (versión rápida, ortonormal) y verifica con una matriz D de 8×8 que $DD^T = I$ y que `dctn` coincide con $DXD^T$; `haar.py` calcula Haar multinivel con PyWavelets y verifica con una matriz H de 8×8 que $HH^T = I$.
2. **Aproximación:** conserva 1, 5, 10 y 25 % de los coeficientes con dos reglas: **soporte prefijado** (esquina de bajas frecuencias / niveles gruesos) y **mayor magnitud**.
3. **Métricas:** error relativo ‖A−Â‖F/‖A‖F, MSE, PSNR y energía retenida; además comprueba que el error cuadrático es igual a la energía descartada.
4. **Cuantización:** q = round(c/Δ) con Δ = 0.05 y costo en bits estimado (cabecera + valores + posiciones). Referencia cruda: 8 bits por píxel.

## Referencias

[1] N. Ahmed, T. Natarajan y K. R. Rao, "Discrete cosine transform," *IEEE Transactions on Computers*, vol. C-23, no. 1, pp. 90–93, ene. 1974, doi: 10.1109/T-C.1974.223784.

[2] G. Strang, "The discrete cosine transform," *SIAM Review*, vol. 41, no. 1, pp. 135–147, 1999, doi: 10.1137/S0036144598336745.

[3] S. G. Mallat, "A theory for multiresolution signal decomposition: The wavelet representation," *IEEE Transactions on Pattern Analysis and Machine Intelligence*, vol. 11, no. 7, pp. 674–693, jul. 1989, doi: 10.1109/34.192463.

[4] O. Velásquez Castañón, "Compresión de imágenes mediante bases ortogonales: Fourier, cosenos, Haar y Walsh–Hadamard; extensión a video," material de referencia, Álgebra Lineal Aplicada 2026-2, Universidad del Pacífico, Lima, oct. 2026.

[5] SciPy, "Discrete Fourier transforms (scipy.fft)," documentación técnica. [En línea]. Disponible: https://docs.scipy.org/doc/scipy/reference/fft.html. [Accedido: 1 nov. 2026].

## Declaración de uso de IA
Herramienta: Claude (Anthropic), modelo Claude Opus 5.5. Finalidad: apoyo en la escritura del código base y la estructura del README. 
