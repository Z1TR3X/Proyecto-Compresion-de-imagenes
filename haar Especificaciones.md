# haar_escala_real.py — Haar por bloques 8×8 con el tamaño real de las imágenes

Versión de Haar análoga a `dct_escala_real.py` del grupo. A diferencia de `haar.py` (256 × 256), aquí **no se reescala**: cada imagen conserva su tamaño real.

| Imagen | Tamaño real (ancho × alto) | Peso original (M·N bytes) |
|---|---|---|
| Resonancia (`Gradiente.webp`) | 728 × 725 | 515,4 KB |
| Radiografía (`bordesmarcados.webp`) | 611 × 626 | 373,5 KB |
| Tomografía (`texturayperiodicidad.jpeg`) | 630 × 425 | 261,5 KB |

## Cómo correrlo
```bash
pip install -r requirements.txt
python haar_escala_real.py
```
Lee `imagenes/` y guarda en `resultados_escala_real/` las figuras (`haar_escala_real_*.png`, subtítulos `% Coef | Error Rel | MSE`) y la tabla (`tabla_haar_escala_real.csv`). No usa PSNR ni cuantización.

## Qué hace (mismo esquema que `dct_escala_real.py`)
1. Carga la imagen en escala de grises de 8 bits (valores 0–255), tamaño real $M \times N$.
2. Padding con `mode="edge"` hasta múltiplos de 8; al final se recorta al tamaño original.
3. Matriz ortonormal explícita **$H_8$** (construida de forma recursiva, $H_N = \operatorname{diag}(H_{N/2}, I_{N/2}) \cdot \text{etapa}(N)$), con verificación $H H^T = I$.
4. Por cada bloque 8×8: $C = H A H^T$; se conservan **$k = \max(1, \lfloor 64 \cdot p \rfloor)$** coeficientes (misma regla que el script de la DCT); $A_{\text{rec}} = H^T M H$. Se verifica Parseval ($\Vert{}A\Vert{}_F^2 = \Vert{}C\Vert{}_F^2$, sumando los bloques) y la reconstrucción exacta al conservar todo.
5. Reglas de selección (por bloque):
   - **prefijada**: esquina superior izquierda del bloque en cuadrados de escala 1, 2, 4, 8 (el último cuadrado puede quedar incompleto). Las posiciones son fijas, no se almacenan.
   - **magnitud**: los $k$ coeficientes de mayor $\vert{}c\vert{}$ de cada bloque; hay que guardar su posición.
6. Métricas: error relativo $\Vert{}A - A_{\text{rec}}\Vert{}_F / \Vert{}A\Vert{}_F$, MSE (en escala 0–255, tras recortar a $[0, 255]$ como en el script de la DCT), % de coeficientes y tamaño estimado (KB) / bpp.

## Decisiones y qué cambia respecto a `haar.py` (256 × 256)
| Aspecto | `haar.py` | `haar_escala_real.py` |
|---|---|---|
| Tamaño | reescalado a 256 × 256 | tamaño real (con padding edge a múltiplos de 8) |
| Transformada | una sola imagen, Haar multinivel | bloques 8 × 8 independientes, $C = H_8 A H_8^T$ ($H_8$ completa, 3 niveles) |
| Selección | $k$ global sobre toda la imagen | **$k$ por bloque** (`magnitud` elige dentro de cada bloque) |
| Escala de píxeles | $[0, 1]$ (MSE pequeño) | 0–255 (MSE ×65 025 respecto a $[0, 1]$; el error relativo no cambia) |
| Tamaño estimado | cabecera + valores con ancho de bits dependiente de los datos + 16 bits por posición | **8 bits por valor** (como `bpp = 8/CR` del script de la DCT) + **6 bits por posición** solo en `magnitud` |

- **Porcentajes:** 1, 5, 10 y 25 %, alineados exactamente con las tasas evaluadas por el grupo.
- **1 % con bloques de 64 coeficientes:** 0,64 coeficientes por bloque, así que, igual que en el script de la DCT, se conserva $\max(1, \lfloor 64 \cdot p \rfloor)$. Por eso los porcentajes reales son: 1 % → $k=1$ (1,56 %), 5 % → $k=3$ (4,69 %), 10 % → $k=6$ (9,38 %), 25 % → $k=16$ (25 %). La tabla muestra `p% (k=…)`. Con $k=1$ ambas reglas conservan prácticamente el coeficiente DC y dan el mismo error (o casi idéntico en Resonancia).
- **Tamaño estimado:** $\text{bloques} \cdot k \cdot (8 [+6 \text{ si magnitud}]) / 8 / 1024$, con los bloques del tamaño con padding. Es una estimación optimista (valores de 8 bits, sin cabecera) pensada para ser comparable con el script de la DCT; no es la misma convención que `haar.py`, así que los KB de ambos scripts no se comparan entre sí.
- **Selección exacta:** se escogen exactamente $k$ coeficientes por bloque (en el script de la DCT un umbral `>=` puede conservar más si hay empates).
- Las dimensiones se informan como ancho × alto (728 × 725); internamente $M = \text{filas} = 725$, $N = \text{columnas} = 728$.

## Resultados
| Imagen | Regla | % Coef | Error Relativo | MSE | Tamaño estimado (KB) | bpp | Peso original (KB) |
|---|---|---|---|---|---|---|---|
| Resonancia (728x725) | prefijada | 1% (k=1) | 0.2553 | 345.45 | 8.1 | 0.126 | 515.4 |
| Resonancia (728x725) | prefijada | 5% (k=3) | 0.1556 | 128.39 | 24.3 | 0.377 | 515.4 |
| Resonancia (728x725) | prefijada | 10% (k=6) | 0.1213 | 78.05 | 48.5 | 0.753 | 515.4 |
| Resonancia (728x725) | prefijada | 25% (k=16) | 0.0698 | 25.86 | 129.4 | 2.008 | 515.4 |
| Resonancia (728x725) | magnitud | 1% (k=1) | 0.2549 | 344.55 | 14.2 | 0.220 | 515.4 |
| Resonancia (728x725) | magnitud | 5% (k=3) | 0.1293 | 88.66 | 42.5 | 0.659 | 515.4 |
| Resonancia (728x725) | magnitud | 10% (k=6) | 0.0810 | 34.81 | 84.9 | 1.318 | 515.4 |
| Resonancia (728x725) | magnitud | 25% (k=16) | 0.0339 | 6.09 | 226.4 | 3.514 | 515.4 |
| Radiografía (611x626) | prefijada | 1% (k=1) | 0.0948 | 191.94 | 5.9 | 0.127 | 373.5 |
| Radiografía (611x626) | prefijada | 5% (k=3) | 0.0727 | 112.81 | 17.8 | 0.382 | 373.5 |
| Radiografía (611x626) | prefijada | 10% (k=6) | 0.0471 | 47.38 | 35.6 | 0.763 | 373.5 |
| Radiografía (611x626) | prefijada | 25% (k=16) | 0.0352 | 26.52 | 95.0 | 2.036 | 373.5 |
| Radiografía (611x626) | magnitud | 1% (k=1) | 0.0948 | 191.94 | 10.4 | 0.223 | 373.5 |
| Radiografía (611x626) | magnitud | 5% (k=3) | 0.0465 | 46.25 | 31.2 | 0.668 | 373.5 |
| Radiografía (611x626) | magnitud | 10% (k=6) | 0.0295 | 18.59 | 62.4 | 1.336 | 373.5 |
| Radiografía (611x626) | magnitud | 25% (k=16) | 0.0145 | 4.50 | 166.3 | 3.562 | 373.5 |
| Tomografía (630x425) | prefijada | 1% (k=1) | 0.1514 | 539.26 | 4.2 | 0.127 | 261.5 |
| Tomografía (630x425) | prefijada | 5% (k=3) | 0.1133 | 301.95 | 12.5 | 0.382 | 261.5 |
| Tomografía (630x425) | prefijada | 10% (k=6) | 0.0948 | 211.39 | 25.0 | 0.765 | 261.5 |
| Tomografía (630x425) | prefijada | 25% (k=16) | 0.0641 | 96.76 | 66.7 | 2.039 | 261.5 |
| Tomografía (630x425) | magnitud | 1% (k=1) | 0.1514 | 539.26 | 7.3 | 0.223 | 261.5 |
| Tomografía (630x425) | magnitud | 5% (k=3) | 0.0936 | 205.84 | 21.9 | 0.669 | 261.5 |
| Tomografía (630x425) | magnitud | 10% (k=6) | 0.0675 | 107.20 | 43.7 | 1.338 | 261.5 |
| Tomografía (630x425) | magnitud | 25% (k=16) | 0.0363 | 31.05 | 116.6 | 3.569 | 261.5 |

Todas las comprobaciones dieron `True`: $H_8 H_8^T = I_8$, Parseval y reconstrucción exacta al 100 %.

## Nota de interpretación
- Con la regla **prefijada** el error baja despacio al subir $k$ (p. ej. Resonancia 0,2553 → 0,0698 de 1 % a 25 %) y deja un cuadriculado por bloque; con **magnitud** baja mucho más rápido (Resonancia 0,2549 → 0,0339; Radiografía 0,0948 → 0,0145; Tomografía 0,1514 → 0,0363), porque cada bloque guarda sus coeficientes más grandes (bordes y texturas).
- El costo es mayor: al 25 % la regla por magnitud pesa 116,6–226,4 KB frente a 66,7–129,4 KB de la prefijada, por los 6 bits de posición por coeficiente.
- A igual tamaño estimado, magnitud al 10 % (Resonancia 84,9 KB, error 0,0810) es comparable a prefijada al 25 % (129,4 KB, error 0,0698): el beneficio de elegir por magnitud se reduce cuando se paga el costo de las posiciones.
- La Tomografía es la imagen más difícil (mayor error con cualquier regla), como en la versión 256 × 256.

