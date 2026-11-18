# Compresión de radiografías con bases ortonormales: DCT vs Haar

**Curso:** Álgebra Lineal Aplicada 2026-2 — Universidad del Pacífico
**Tema:** Compresión de imágenes mediante bases ortogonales
**Video:** [Ver video](https://drive.google.com/drive/folders/1949UsEJW-8gNVBSwZbvh5TFSar30c8Aw?usp=sharing)

## Integrantes
- Kevin Paolo Julca Lucano
- Tania Lisset Chavez Alvitez
- Jose Julian Tucto Vargas
- Rosa Leticia Ortiz Palomino


## Problema y objetivo
Comparar qué base ortonormal (DCT global o Haar multinivel) conserva mejor la información de distintos tipos de imágenes médicas (suave, con bordes y con textura) cuando se guarda solo una fracción de los coeficientes, medido con el error de recuperación, el MSE, la tasa de compresión y el tamaño comprimido estimado en bits.


## Organización de archivos
```
dct.py                -> experimento con la DCT global (tamaño real de la imagen)
haar.py               -> experimento con Haar multinivel
Dct.md                -> desarrollo matemático de la DCT-II
Haar (1).md           -> desarrollo matemático de Haar
requirements.txt      -> dependencias
imagenes/             -> imágenes de entrada
resultados/           -> figuras e imágenes comprimidas que generan los scripts
```


## Cómo ejecutarlo
Requisitos: Python 3.10 o superior.

```bash
pip install -r requirements.txt
python dct.py                   # experimento con DCT
python haar.py                  # experimento con Haar
```

### Opciones de `dct.py`
Se cambian directamente al inicio del archivo:

| Opción | Qué hace |
|---|---|
| `RUTA` | Imagen a usar, por ejemplo `"imagenes/Gradiente.webp"`. Se procesa una imagen por ejecución. |
| `PORCENTAJES` | Porcentajes de coeficientes que se conservan: `[1, 5, 10, 25]`, como pide el experimento. |
| `DELTA` | Paso de cuantización (Δ = 0.05). |
| `GUARDAR_PREFIJADA` | `True` = guarda las imágenes comprimidas con la regla prefijada. |
| `GUARDAR_MAGNITUD` | `True` = guarda las imágenes comprimidas con la regla por magnitud. |
| `GUARDAR_CUADRO` | `True` = guarda el cuadro comparativo (prefijada vs magnitud). |
| `MOSTRAR_SIN_DELTA` | `True` = agrega a la tabla el error sin Δ, el MSE sin Δ y la energía retenida. |

La tabla de resultados siempre se imprime en la terminal.

`dct.py` trabaja con la imagen **en su tamaño real**: solo la convierte a escala de grises y la normaliza al rango [0, 1], sin recortarla ni reescalarla. `haar.py` todavía recorta un cuadrado central y lo lleva a 256×256.

## Qué hacen `dct.py` y `haar.py`
Los dos archivos siguen los mismos pasos, cada uno con su base: la DCT [1], [2] y Haar ortonormal multinivel [3]. El procedimiento sigue el material del curso [4].

1. **Verificaciones:** `dct.py` calcula la DCT con `scipy.fft.dctn` [5] (versión rápida, ortonormal) y verifica con una matriz D de 8×8 que $DD^T = I$ y que `dctn` coincide con $DXD^T$; además comprueba Parseval ($\|A\|^2 = \|C\|^2$) en la imagen. `haar.py` calcula Haar multinivel con PyWavelets, verifica con una matriz H de 8×8 que $HH^T = I$ y comprueba la reconstrucción exacta y Parseval.
2. **Aproximación:** conserva 1, 5, 10 y 25 % de los coeficientes con dos reglas: **soporte prefijado** (esquina de bajas frecuencias / niveles gruesos) y **mayor magnitud**.
3. **Cuantización:** q = round(c/Δ) con Δ = 0.05, y se reconstruye a partir de Δ·q.
4. **Métricas:**
   - **MSE:** desviación promedio por píxel, $\frac{1}{mn}\sum (A-\hat A)^2$.
   - **Error de recuperación:** $\|A-\hat A\|_F / \|A\|_F$, en %.
   - **Tasa de compresión:** bits de la imagen cruda / bits comprimidos.
   - **Tamaño comprimido:** bits estimados = cabecera (64) + valores + posiciones (solo con la regla por magnitud).
   - Referencia cruda: 8 bits por píxel. No se compara contra el peso de los archivos JPG/WEBP de entrada, que ya vienen comprimidos.

## Resultados principales (DCT, tamaño real)
Error de recuperación y tasa de compresión al conservar el 10 % de los coeficientes:

| Imagen | Tamaño real | Prefijada | Magnitud |
|---|---|---|---|
| Gradiente (suave) | 725 × 728 | 3.83 % · 6.2:1 | 3.04 % · 2.4:1 |
| Bordesmarcados (bordes) | 626 × 611 | 3.85 % · 5.7:1 | 2.38 % · 2.4:1 |
| Texturayperiodicidad (textura) | 425 × 630 | 6.65 % · 5.7:1 | 5.01 % · 2.4:1 |

- Con el mismo porcentaje, la regla **por magnitud siempre tiene menor error**.
- Pero **menor error no significa menos bits**: por magnitud hay que guardar la posición de cada coeficiente. Al 25 %, la regla por magnitud tiene una tasa de 1.0:1 (no comprime), mientras que la prefijada todavía comprime entre 2.3:1 y 2.5:1.
- La imagen con **textura** es la más difícil de comprimir: los detalles finos reparten la energía en muchos coeficientes.
- Limitación: son solo tres imágenes, así que no se concluye que una base o regla sea siempre mejor.

## Comparación DCT vs Haar (tamaño real)
Se comparó qué base pierde menos información cuando la imagen comprimida pesa lo mismo (DCT sobre toda la imagen; Haar por bloques de 8 × 8):

| Imagen | Tipo | Mejor base |
|---|---|---|
| Resonancia | Gradiente suave | **DCT**: gana en todos los tamaños (≈ 84 KB: 3,83 % frente a 8,10 % de error) |
| Radiografía | Bordes marcados | **Haar**: gana desde ≈ 30 KB (≈ 64 KB: 2,95 % frente a 3,85 %) |
| Tomografía | Textura | **Empate**: con poco espacio gana la DCT; con más espacio, Haar |

La mejor base **sí depende del tipo de imagen**. Limitación: DCT y Haar no usan la misma división de la imagen (completa vs bloques de 8 × 8).

## Referencias

[1] N. Ahmed, T. Natarajan y K. R. Rao, "Discrete cosine transform," *IEEE Transactions on Computers*, vol. C-23, no. 1, pp. 90–93, ene. 1974, doi: 10.1109/T-C.1974.223784.

[2] G. Strang, "The discrete cosine transform," *SIAM Review*, vol. 41, no. 1, pp. 135–147, 1999, doi: 10.1137/S0036144598336745.

[3] S. G. Mallat, "A theory for multiresolution signal decomposition: The wavelet representation," *IEEE Transactions on Pattern Analysis and Machine Intelligence*, vol. 11, no. 7, pp. 674–693, jul. 1989, doi: 10.1109/34.192463.

[4] SciPy, "Discrete Fourier transforms (scipy.fft)," documentación técnica. [En línea]. Disponible: https://docs.scipy.org/doc/scipy/reference/fft.html. [Accedido: 1 nov. 2026].



## Declaración de uso de IA
- **Herramienta:** Claude (Anthropic)
- **Finalidad y partes del trabajo:** apoyo en la escritura del código base de `dct.py` y `haar.py`,
