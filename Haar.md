# Desarrollo matemático de la Transformada de Haar

**Proyecto:** Comparación de la DCT-II y la Transformada de Haar ortonormal para la representación de imágenes en escala de grises.

## 1. Fundamento matemático de la Transformada de Haar

La Transformada de Haar es una transformada wavelet ortonormal que representa una señal mediante coeficientes de aproximación y de detalle. Su construcción utiliza funciones constantes por tramos y combinaciones locales de suma y diferencia. La aproximación describe el comportamiento general de una señal en una región, mientras que los detalles representan las variaciones que aparecen dentro de ella [1], [2].

La transformada se fundamenta en dos funciones: la **función de escala**, asociada a la aproximación, y la **wavelet madre**, asociada al detalle. A partir de estas funciones se construyen bases ortonormales que permiten analizar una señal en diferentes escalas.

### 1.1. Función de escala de Haar

La función de escala de Haar, denotada por `φ(x)`, se define como:

$$
\phi(x)=
\begin{cases}
1, & 0\leq x<1,\\
0, & \text{en otro caso}.
\end{cases}
$$

La función es constante en el intervalo `[0,1)` y nula fuera de él. Representa la componente de aproximación de la señal dentro del intervalo considerado.

### 1.2. Wavelet madre de Haar

La wavelet madre, denotada por `ψ(x)`, se define como:

$$
\psi(x)=
\begin{cases}
1, & 0\leq x<\frac{1}{2},\\
-1, & \frac{1}{2}\leq x<1,\\
0, & \text{en otro caso}.
\end{cases}
$$

La wavelet presenta un tramo positivo y otro negativo. Esta estructura permite medir la diferencia entre la primera y la segunda mitad del intervalo, por lo que representa la información de detalle.

Ambas funciones tienen norma unitaria:

$$
\int_{-\infty}^{\infty}\phi(x)^2\,dx=1,
\qquad
\int_{-\infty}^{\infty}\psi(x)^2\,dx=1.
$$

Además, son ortogonales:

$$
\int_{-\infty}^{\infty}\phi(x)\psi(x)\,dx=0.
$$

La integral es cero porque las dos mitades de la wavelet tienen la misma longitud y signos opuestos. Estas propiedades son la base de la construcción ortonormal de Haar [1], [3].

## 2. Construcción de la Transformada de Haar unidimensional

### 2.1. Transformación de dos valores

Sean `x₁` y `x₂` dos componentes de una señal. La transformación de Haar produce un coeficiente de aproximación `a` y uno de detalle `d`:

$$
a=\frac{x_1+x_2}{\sqrt{2}},
\qquad
d=\frac{x_1-x_2}{\sqrt{2}}.
$$

En forma matricial:

$$
\begin{pmatrix}
a\\
d
\end{pmatrix}
=
H_2
\begin{pmatrix}
x_1\\
x_2
\end{pmatrix},
\qquad
H_2=
\frac{1}{\sqrt{2}}
\begin{pmatrix}
1&1\\
1&-1
\end{pmatrix}.
$$

La primera fila de `H₂` calcula la suma normalizada de los dos valores; la segunda calcula su diferencia normalizada. El factor `1/√2` permite que ambas filas tengan norma unitaria.

Si `x=[x₁,x₂]ᵀ` es el vector original y `c=[a,d]ᵀ` el vector transformado, entonces la transformación se expresa como `c=H₂x`. Es una transformación lineal e invertible.

### 2.2. Demostración de ortonormalidad de `H₂`

Una matriz real cuadrada `H` es ortogonal si el producto de la matriz por su transpuesta es la identidad:

$$
HH^T=I.
$$

Para la matriz de Haar de orden dos:

$$
\begin{aligned}
H_2H_2^T
&=
\frac{1}{2}
\begin{pmatrix}
1&1\\
1&-1
\end{pmatrix}
\begin{pmatrix}
1&1\\
1&-1
\end{pmatrix}\\
&=
\frac{1}{2}
\begin{pmatrix}
2&0\\
0&2
\end{pmatrix}
=
\begin{pmatrix}
1&0\\
0&1
\end{pmatrix}
=I_2.
\end{aligned}
$$

Por tanto, las filas de `H₂` son ortogonales entre sí y tienen norma uno. Se concluye que:

$$
H_2^{-1}=H_2^T.
$$

Como `H₂` es simétrica, también se cumple `H₂ᵀ=H₂`. La transformación inversa es:

$$
\begin{pmatrix}
x_1\\
x_2
\end{pmatrix}
=
H_2^T
\begin{pmatrix}
a\\
d
\end{pmatrix}.
$$

### 2.3. Construcción de la matriz `H₄`

Para cuatro componentes se aplica el procedimiento de manera jerárquica. Primero se forman dos pares y se calculan sus aproximaciones y detalles. Después, las dos aproximaciones se vuelven a transformar para obtener una aproximación global y un detalle de escala intermedia.

Con el orden de coeficientes: aproximación global, detalle de escala intermedia y dos detalles locales, se obtiene:

$$
H_4=
\begin{pmatrix}
\frac{1}{2}&\frac{1}{2}&\frac{1}{2}&\frac{1}{2}\\
\frac{1}{2}&\frac{1}{2}&-\frac{1}{2}&-\frac{1}{2}\\
\frac{1}{\sqrt{2}}&-\frac{1}{\sqrt{2}}&0&0\\
0&0&\frac{1}{\sqrt{2}}&-\frac{1}{\sqrt{2}}
\end{pmatrix}.
$$

La primera fila representa la aproximación conjunta de los cuatro componentes. La segunda compara la primera mitad con la segunda. Las dos últimas filas representan las diferencias dentro de cada par.

La matriz es ortonormal:

$$
H_4H_4^T=I_4,
\qquad
H_4^{-1}=H_4^T.
$$

Este procedimiento muestra que los coeficientes de detalle se generan en distintas escalas, mientras que la aproximación continúa descomponiéndose en el siguiente nivel.

### 2.4. Generalización a una señal de longitud `N=2^J`

El procedimiento se extiende a señales cuya longitud es una potencia de dos, `N=2^J`, donde `J` es un entero no negativo. En cada nivel se agrupan los coeficientes de aproximación de dos en dos. Las ecuaciones de análisis son:

$$
a_{j+1,k}
=
\frac{a_{j,2k}+a_{j,2k+1}}{\sqrt{2}},
$$

$$
d_{j+1,k}
=
\frac{a_{j,2k}-a_{j,2k+1}}{\sqrt{2}}.
$$

Aquí, `a_j,k` representa un coeficiente de aproximación en el nivel `j`, mientras que `d_(j+1),k` representa el detalle generado en el siguiente nivel. El índice `k` identifica cada pareja.

En cada iteración, la transformación se aplica únicamente a los coeficientes de aproximación. Los detalles de los niveles anteriores se conservan. Al finalizar, la señal queda representada por una aproximación global y por detalles asociados a diferentes escalas [1], [3].

### 2.5. Estructura general de la matriz de Haar `H_N`

La matriz ortonormal de Haar de orden `N=2^J` se construye con una fila de aproximación global y filas de detalle con diferentes longitudes de soporte.

La fila de aproximación global es:

$$
h_0=\frac{1}{\sqrt{N}}
\begin{pmatrix}
1&1&\cdots&1
\end{pmatrix}.
$$

Para describir una fila de detalle, sea `L=2^r` su longitud de soporte. La fila contiene `L/2` valores positivos seguidos por `L/2` valores negativos, todos con magnitud `1/√L`. Las demás posiciones son cero. Formalmente:

$$
h_{r,k}[n]=
\begin{cases}
\frac{1}{\sqrt{L}}, & kL\leq n<kL+\frac{L}{2},\\
-\frac{1}{\sqrt{L}}, & kL+\frac{L}{2}\leq n<(k+1)L,\\
0, & \text{en otro caso},
\end{cases}
$$

donde `L=2^r`, `k` identifica el bloque y `n` la posición de la señal. Para cada longitud `L`, los bloques no se superponen y cubren la señal.

La matriz `H_N` se forma ordenando estas filas por escalas. Una convención frecuente coloca primero la aproximación global, después los detalles de mayor soporte y finalmente los detalles locales. El orden puede variar entre implementaciones; una permutación de filas cambia la ubicación de los coeficientes, pero conserva la ortonormalidad.

Para `N=8`, la matriz correspondiente a esta convención es:

$$
H_8=
\begin{pmatrix}
\frac{1}{\sqrt{8}}&\frac{1}{\sqrt{8}}&\frac{1}{\sqrt{8}}&\frac{1}{\sqrt{8}}&\frac{1}{\sqrt{8}}&\frac{1}{\sqrt{8}}&\frac{1}{\sqrt{8}}&\frac{1}{\sqrt{8}}\\
\frac{1}{\sqrt{8}}&\frac{1}{\sqrt{8}}&\frac{1}{\sqrt{8}}&\frac{1}{\sqrt{8}}&-\frac{1}{\sqrt{8}}&-\frac{1}{\sqrt{8}}&-\frac{1}{\sqrt{8}}&-\frac{1}{\sqrt{8}}\\
\frac{1}{2}&\frac{1}{2}&-\frac{1}{2}&-\frac{1}{2}&0&0&0&0\\
0&0&0&0&\frac{1}{2}&\frac{1}{2}&-\frac{1}{2}&-\frac{1}{2}\\
\frac{1}{\sqrt{2}}&-\frac{1}{\sqrt{2}}&0&0&0&0&0&0\\
0&0&\frac{1}{\sqrt{2}}&-\frac{1}{\sqrt{2}}&0&0&0&0\\
0&0&0&0&\frac{1}{\sqrt{2}}&-\frac{1}{\sqrt{2}}&0&0\\
0&0&0&0&0&0&\frac{1}{\sqrt{2}}&-\frac{1}{\sqrt{2}}
\end{pmatrix}.
$$

Las dos primeras filas representan la aproximación global y el contraste entre las dos mitades. Las filas tercera y cuarta representan diferencias entre grupos de cuatro componentes. Las últimas cuatro filas representan diferencias locales entre pares consecutivos.

## 3. Propiedades de ortonormalidad

La ortonormalidad significa que cada fila de la matriz tiene norma uno y que el producto interno entre dos filas distintas es cero. Para una matriz de Haar cuadrada `H_N`, ambas condiciones se expresan mediante:

$$
H_NH_N^T=I_N.
$$

### 3.1. Norma unitaria de las filas

La norma euclidiana al cuadrado de una fila `h_i` se define como:

$$
\|h_i\|^2=\sum_{n=0}^{N-1}h_i[n]^2.
$$

La fila de aproximación global tiene `N` componentes iguales a `1/√N`. Por tanto:

$$
\|h_0\|^2
=
N\left(\frac{1}{\sqrt{N}}\right)^2
=
N\left(\frac{1}{N}\right)=1.
$$

Una fila de detalle con soporte de longitud `L` tiene `L` componentes no nulos, cada uno con magnitud `1/√L`. Entonces:

$$
\|h_{r,k}\|^2
=
L\left(\frac{1}{\sqrt{L}}\right)^2
=
L\left(\frac{1}{L}\right)=1.
$$

Así, todas las filas tienen norma unitaria.

### 3.2. Ortogonalidad entre filas distintas

El producto interno entre dos filas `h_p` y `h_q` se define como:

$$
\langle h_p,h_q\rangle
=
\sum_{n=0}^{N-1}h_p[n]h_q[n].
$$

En la base de Haar, dos filas de detalle tienen soportes disjuntos o uno de sus soportes está contenido en una mitad constante del otro. En el primer caso, el producto interno es cero porque no comparten posiciones no nulas. En el segundo, los valores positivos y negativos de la fila de menor soporte se cancelan al sumarse.

La fila de aproximación global también es ortogonal a todas las filas de detalle, ya que cada fila de detalle contiene la misma cantidad de valores positivos y negativos.

Por ejemplo, el producto interno entre las dos primeras filas de `H₈` es:

$$
\begin{aligned}
\langle h_0,h_1\rangle
&=
4\left(\frac{1}{\sqrt{8}}\frac{1}{\sqrt{8}}\right)
+
4\left(\frac{1}{\sqrt{8}}\frac{-1}{\sqrt{8}}\right)\\
&=\frac{4}{8}-\frac{4}{8}=0.
\end{aligned}
$$

Por tanto, las filas son mutuamente ortogonales. Al combinar esta propiedad con la norma unitaria se obtiene:

$$
\boxed{H_NH_N^T=I_N}.
$$

En consecuencia, la inversa de la matriz de Haar es su transpuesta:

$$
H_N^{-1}=H_N^T.
$$

## 4. Transformada de Haar bidimensional

Una imagen digital en escala de grises se representa mediante una matriz:

$$
A\in\mathbb{R}^{M\times N},
$$

donde `a_ij` representa la intensidad del píxel ubicado en la fila `i` y la columna `j`.

La extensión de Haar a dos dimensiones se obtiene aplicando la transformación unidimensional sobre ambos ejes. Esta propiedad se denomina **separabilidad**: la transformación bidimensional se calcula mediante dos productos matriciales consecutivos [3], [4].

Sean `H_M` la matriz de Haar de orden `M` y `H_N` la matriz de Haar de orden `N`. La transformación directa es:

$$
C=H_MAH_N^T.
$$

Puede calcularse en dos etapas:

$$
B=H_MA,
$$

$$
C=BH_N^T.
$$

La primera multiplicación transforma un eje de la imagen y la segunda transforma el otro. La matriz `C` contiene los coeficientes de aproximación y detalle en ambas dimensiones.

Si la imagen es cuadrada, con `M=N`, se utiliza la misma matriz en ambos ejes:

$$
C=H_NAH_N^T.
$$

### 4.1. Transformación inversa bidimensional

Por la ortonormalidad de las matrices se cumple:

$$
H_M^{-1}=H_M^T,
\qquad
H_N^{-1}=H_N^T.
$$

Partiendo de la transformación directa:

$$
C=H_MAH_N^T,
$$

se multiplica por `H_Mᵀ` a la izquierda y por `H_N` a la derecha:

$$
H_M^TCH_N
=
H_M^TH_MAH_N^TH_N.
$$

Como:

$$
H_M^TH_M=I_M,
\qquad
H_N^TH_N=I_N,
$$

se obtiene la reconstrucción:

$$
\boxed{A=H_M^TCH_N}.
$$

La igualdad demuestra que la transformación bidimensional es reversible cuando se conservan todos los coeficientes.

## 5. Subbandas de la transformada bidimensional

En una etapa de la transformación bidimensional, cada eje se separa en componentes de aproximación (L) y de detalle (H). Al combinar ambas posibilidades se obtienen cuatro subbandas:

- **LL:** aproximación en ambos ejes. Representa la componente general o suavizada de la imagen.
- **LH:** aproximación en un eje y detalle en el otro.
- **HL:** detalle en un eje y aproximación en el otro.
- **HH:** detalle en ambos ejes.

Matemáticamente, estas subbandas se relacionan con los productos de la función de escala `φ` y la wavelet madre `ψ`:

$$
\phi(x)\phi(y),\qquad
\phi(x)\psi(y),\qquad
\psi(x)\phi(y),\qquad
\psi(x)\psi(y).
$$

La primera combinación representa aproximación en las dos coordenadas. Las otras tres representan variaciones en uno o ambos ejes. La interpretación horizontal o vertical de LH y HL depende de qué eje se transforma primero y de la convención utilizada por la implementación [4].

En una descomposición multinivel, el procedimiento se aplica nuevamente sobre la subbanda LL. De esta manera, la aproximación de un nivel se descompone en una nueva aproximación y tres grupos adicionales de detalles. Se obtiene así una representación jerárquica de la información [1], [3].

## 6. Análisis multirresolución

El análisis multirresolución permite interpretar la descomposición wavelet como una secuencia de espacios de aproximación anidados:

$$
\cdots\subset V_{-1}\subset V_0\subset V_1\subset\cdots
$$

Cada espacio `V_j` contiene aproximaciones de una función a una resolución determinada. El espacio `W_j` representa la información adicional necesaria para pasar de una resolución a la siguiente. En el análisis ortogonal se cumple:

$$
V_{j+1}=V_j\oplus W_j.
$$

El símbolo `⊕` representa una suma directa ortogonal: los elementos de `V_j` y `W_j` son ortogonales entre sí y juntos generan `V_(j+1)`.

En la construcción de Haar, la función de escala genera los espacios de aproximación y la wavelet madre genera los espacios de detalle. Sus versiones dilatadas y trasladadas se expresan como:

$$
\phi_{j,k}(x)=2^{j/2}\phi(2^jx-k),
$$

$$
\psi_{j,k}(x)=2^{j/2}\psi(2^jx-k).
$$

Aquí, `j` identifica la escala y `k` la traslación. El factor `2^(j/2)` normaliza las funciones para conservar su norma unitaria. Este marco explica cómo una señal puede representarse mediante una aproximación en una escala y detalles correspondientes a otras escalas [1], [2].

## 7. Conservación de energía e identidad de Parseval

Una transformación ortonormal conserva la norma euclidiana. Si `x` es el vector original y `c=H_Nx` su representación en la base de Haar, entonces:

$$
\begin{aligned}
\|c\|_2^2
&=c^Tc\\
&=(H_Nx)^T(H_Nx)\\
&=x^TH_N^TH_Nx\\
&=x^Tx\\
&=\|x\|_2^2.
\end{aligned}
$$

Por tanto:

$$
\boxed{\|c\|_2^2=\|x\|_2^2}.
$$

Esta igualdad corresponde a la identidad de Parseval: la suma de los cuadrados de los valores originales es igual a la suma de los cuadrados de los coeficientes transformados.

Para matrices, la misma propiedad se expresa mediante la norma de Frobenius. Si `C=H_MAH_Nᵀ`, entonces:

$$
\boxed{\|C\|_F^2=\|A\|_F^2}.
$$

La transformación redistribuye la energía entre los coeficientes, pero no modifica su cantidad total mientras se conserven todos ellos [1], [3].

## 8. Reconstrucción y relación del error con los coeficientes

La transformación completa es reversible. Sin embargo, si algunos coeficientes se modifican o se sustituyen por cero, la reconstrucción puede diferir de la matriz original.

Sea `C` la matriz completa de coeficientes y `C̃` una matriz modificada. La reconstrucción correspondiente es:

$$
\widetilde{A}=H_M^T\widetilde{C}H_N.
$$

La diferencia entre la matriz original y la reconstruida es:

$$
A-\widetilde{A}
=
H_M^T(C-\widetilde{C})H_N.
$$

Debido a que las matrices de transformación son ortonormales, la norma de Frobenius se conserva. Por tanto:

$$
\boxed{
\|A-\widetilde{A}\|_F^2
=
\|C-\widetilde{C}\|_F^2
}.
$$

Esta igualdad establece que el error cuadrático en el dominio original es igual al error cuadrático entre los coeficientes originales y los modificados en el dominio transformado. Si los coeficientes descartados se reemplazan por cero, el error corresponde a la suma de sus cuadrados.

Esta propiedad fundamenta matemáticamente el análisis de selección de coeficientes. No obstante, la transformada de Haar por sí sola no define un archivo comprimido: la reducción efectiva del almacenamiento también depende de la cuantización, la codificación y la información necesaria para identificar los coeficientes conservados.

## 9. Referencias bibliográficas

[1] S. G. Mallat, “Multiresolution Approximations and Wavelet Orthonormal Bases of `L²(ℝ)`,” *Transactions of the American Mathematical Society*, vol. 315, no. 1, pp. 69–87, 1989, doi: 10.1090/S0002-9947-1989-1008470-5.

[2] A. Haar, “Zur Theorie der orthogonalen Funktionensysteme,” *Mathematische Annalen*, vol. 69, pp. 331–371, 1910, doi: 10.1007/BF01456326.

[3] S. Mallat, “A Theory for Multiresolution Signal Decomposition: The Wavelet Representation,” *IEEE Transactions on Pattern Analysis and Machine Intelligence*, vol. 11, no. 7, pp. 674–693, 1989, doi: 10.1109/34.192463.

[4] I. Daubechies, *Ten Lectures on Wavelets*. Philadelphia, PA, USA: Society for Industrial and Applied Mathematics, 1992, doi: 10.1137/1.9781611970104.
