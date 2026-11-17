# Desarrollo matemático de la Transformada de Haar

**Proyecto:** Comparación de la DCT-II y la Transformada de Haar
ortonormal para la representación de imágenes en escala de grises.

## 1. Fundamento matemático de la Transformada de Haar

La Transformada de Haar es una transformada wavelet ortonormal que
representa una señal mediante coeficientes de aproximación y de detalle.
Su construcción se basa en funciones constantes por tramos y en
operaciones locales de suma y diferencia. A diferencia de una
representación que describe únicamente los valores originales, la base
de Haar permite expresar la señal según su comportamiento en distintas
escalas \[1\], \[2\].

La idea matemática central consiste en separar dos tipos de información:

-   **Aproximación:** describe el comportamiento general de un conjunto
    de datos.
-   **Detalle:** representa las diferencias o variaciones entre partes
    de ese conjunto.

La transformación se construye mediante funciones de escala y funciones
wavelet. Estas funciones forman una familia ortonormal y permiten
descomponer una señal en componentes que pueden analizarse por niveles
de resolución \[1\].

### 1.1. Función de escala de Haar

La función de escala de Haar, denotada por $\phi(x)$, se define
como:

```math
\phi(x)= \begin{cases} 1, & 0\leq x<1,\\ 0, & \text{en otro caso}. \end{cases}
```

Esta función es constante en el intervalo $[0,1)$ y nula fuera de él.
Su función representa el componente de aproximación, porque conserva el
comportamiento promedio de la señal en el intervalo considerado.

### 1.2. Wavelet madre de Haar

La wavelet madre, denotada por $\psi(x)$, se define mediante:

```math
\psi(x)= \begin{cases} 1, & 0\leq x<\frac{1}{2},\\ -1, & \frac{1}{2}\leq x<1,\\ 0, & \text{en otro caso}. \end{cases}
```

La función presenta dos tramos: uno positivo y otro negativo. Esta
diferencia de signos permite medir el contraste entre la primera y la
segunda mitad del intervalo. Por ello, la wavelet madre representa la
información de detalle.

Ambas funciones están normalizadas y son ortogonales. En particular:

```math
\int_{-\infty}^{\infty}\phi(x)^2\,dx=1, \qquad \int_{-\infty}^{\infty}\psi(x)^2\,dx=1,
```

y

```math
\int_{-\infty}^{\infty}\phi(x)\psi(x)\,dx=0.
```

La última igualdad se cumple porque las dos mitades de $\psi(x)$
tienen la misma longitud y signos opuestos. Estas propiedades
constituyen la base de la ortonormalidad de Haar \[1\], \[3\].

## 2. Construcción de la Transformada de Haar unidimensional

### 2.1. Transformación de dos valores

Sean $x_1$ y $x_2$ dos componentes de una señal. La transformación de
Haar los convierte en un coeficiente de aproximación $a$ y uno de
detalle $d$:

```math
a=\frac{x_1+x_2}{\sqrt{2}}, \qquad d=\frac{x_1-x_2}{\sqrt{2}}.
```

Estas expresiones pueden escribirse como un producto matricial:

```math
\begin{bmatrix} a\\ d \end{bmatrix} = H_2 \begin{bmatrix} x_1\\ x_2 \end{bmatrix}, \qquad H_2= \frac{1}{\sqrt{2}} \begin{bmatrix} 1&1\\ 1&-1 \end{bmatrix}.
```

La primera fila de $H_2$ combina los dos valores mediante una suma
normalizada; la segunda los combina mediante una diferencia normalizada.
El factor $1/\sqrt{2}$ garantiza que ambas filas tengan norma
unitaria.

La transformación es lineal: si $x$ es el vector original y $c$ el
vector transformado, entonces $c=H_2x$. Por ser una matriz cuadrada e
invertible, permite recuperar los datos originales cuando se conservan
ambos coeficientes.

### 2.2. Ortonormalidad de la matriz $H_2$

Una matriz cuadrada real $H$ es ortogonal si su producto por la
transpuesta es la matriz identidad:

```math
HH^T=I.
```

Para $H_2$, se calcula:

```math
\begin{aligned} H_2H_2^T &= \frac{1}{2} \begin{bmatrix} 1&1\\ 1&-1 \end{bmatrix} \begin{bmatrix} 1&1\\ 1&-1 \end{bmatrix}\\ &= \frac{1}{2} \begin{bmatrix} 2&0\\ 0&2 \end{bmatrix} = \begin{bmatrix} 1&0\\ 0&1 \end{bmatrix} =I_2. \end{aligned}
```

Por tanto, las filas de $H_2$ son ortogonales entre sí y tienen norma
igual a uno. De esta propiedad se obtiene:

```math
H_2^{-1}=H_2^T.
```

Además, como $H_2$ es simétrica, $H_2^T=H_2$. La transformación inversa
queda expresada como:

```math
\begin{bmatrix} x_1\\ x_2 \end{bmatrix} = H_2^T \begin{bmatrix} a\\ d \end{bmatrix}.
```

### 2.3. Construcción de la matriz $H_4$

La transformación puede extenderse a cuatro componentes aplicando el
mismo procedimiento por niveles. Primero se forman dos pares y se
calculan sus aproximaciones y detalles. Luego, las dos aproximaciones se
vuelven a transformar para obtener una aproximación global y un detalle
de escala intermedia.

Con el orden de coeficientes: aproximación global, detalle de escala
intermedia y dos detalles locales, la matriz resultante es:

```math
H_4= \begin{bmatrix} \frac{1}{2}&\frac{1}{2}&\frac{1}{2}&\frac{1}{2}\\ \frac{1}{2}&\frac{1}{2}&-\frac{1}{2}&-\frac{1}{2}\\ \frac{1}{\sqrt{2}}&-\frac{1}{\sqrt{2}}&0&0\\ 0&0&\frac{1}{\sqrt{2}}&-\frac{1}{\sqrt{2}} \end{bmatrix}.
```

La primera fila representa la aproximación conjunta de los cuatro
componentes. La segunda compara las dos mitades del vector. Las dos
últimas filas registran las diferencias dentro de cada par.

La matriz mantiene la propiedad de ortonormalidad:

```math
H_4H_4^T=I_4.
```

Por consiguiente, su inversa también coincide con su transpuesta:

```math
H_4^{-1}=H_4^T.
```

La construcción muestra que los coeficientes de detalle se generan en
diferentes escalas, mientras que la aproximación se vuelve a descomponer
en el siguiente nivel.

### 2.4. Generalización a una dimensión $N=2^J$

El procedimiento se puede extender a una señal de longitud $N=2^J$,
donde $J$ es un entero no negativo. En cada nivel se agrupan los
coeficientes de aproximación de dos en dos. Para el nivel $j$, las
relaciones de análisis son:

```math
a_{j+1,k} = \frac{a_{j,2k}+a_{j,2k+1}}{\sqrt{2}},
```

```math
d_{j+1,k} = \frac{a_{j,2k}-a_{j,2k+1}}{\sqrt{2}}.
```

Aquí, $a_{j,k}$ representa los coeficientes de aproximación del nivel
$j$, mientras que $d_{j+1,k}$ representa los coeficientes de detalle
generados en el siguiente nivel. El índice $k$ identifica cada pareja de
coeficientes.

La operación se repite únicamente sobre los coeficientes de
aproximación. Los detalles obtenidos en cada nivel se conservan como
información de las variaciones correspondientes a esa escala. De esta
manera, una señal queda representada por una aproximación final y por un
conjunto de detalles de distintas resoluciones \[1\], \[3\].

### 2.5. Estructura general de las filas de $H_N$

Para $N=2^J$, la matriz ortonormal de Haar puede organizarse mediante
una fila de aproximación global y filas de detalle con diferentes
soportes.

La fila de aproximación global es:

```math
h_0= \frac{1}{\sqrt{N}} \begin{bmatrix} 1&1&\cdots&1 \end{bmatrix}.
```

Una fila de detalle con soporte de longitud $L=2^r$ contiene $L/2$
valores positivos, seguidos por $L/2$ valores negativos, todos de
magnitud $1/\sqrt{L}$; las demás posiciones son cero. Los
soportes se desplazan en bloques no superpuestos de longitud $L$.

En forma esquemática, cada fila de detalle tiene la estructura:

```math
h_{r,k}[n]= \begin{cases} \frac{1}{\sqrt{L}}, & kL\leq n<kL+\frac{L}{2},\\ -\frac{1}{\sqrt{L}}, & kL+\frac{L}{2}\leq n<(k+1)L,\\ 0, & \text{en otro caso}, \end{cases}
```

donde $L=2^r$, $k$ identifica el bloque y $n$ identifica la posición de
la señal. La fila se normaliza con $1/\sqrt{L}$, de modo que su
norma sea uno.

La matriz $H_N$ se forma colocando estas filas en un orden por escalas.
Una convención frecuente ubica primero la aproximación global, luego los
detalles de mayor soporte y finalmente los detalles locales. El orden de
las filas puede variar entre textos o implementaciones; una permutación
de filas cambia el orden de los coeficientes, pero conserva la
ortonormalidad.

## 3. Ortonormalidad de la matriz de Haar

La ortonormalidad significa que cada vector de la base tiene norma uno y
que dos vectores diferentes tienen producto interno igual a cero. Para
una matriz de Haar cuadrada $H_N$, estas condiciones se expresan
conjuntamente como:

```math
H_NH_N^T=I_N.
```

### 3.1. Norma unitaria de las filas

La norma euclidiana al cuadrado de una fila $h_i$ se define como:

```math
\|h_i\|^2=\sum_{n=0}^{N-1}h_i[n]^2.
```

Para la fila de aproximación global, todos sus $N$ componentes son
$1/\sqrt{N}$. Por tanto:

```math
\|h_0\|^2 = N\left(\frac{1}{\sqrt{N}}\right)^2 = N\left(\frac{1}{N}\right)=1.
```

Para una fila de detalle con soporte de longitud $L$, existen $L$
componentes no nulos, cada uno con magnitud $1/\sqrt{L}$.
Entonces:

```math
\|h_{r,k}\|^2 = L\left(\frac{1}{\sqrt{L}}\right)^2 = L\left(\frac{1}{L}\right)=1.
```

Así, todas las filas de la matriz tienen norma unitaria.

### 3.2. Ortogonalidad entre filas distintas

El producto interno entre dos filas $h_p$ y $h_q$ se define como:

```math
\langle h_p,h_q\rangle = \sum_{n=0}^{N-1}h_p[n]h_q[n].
```

Para las filas de detalle de Haar se presentan dos situaciones:

1.  **Soportes disjuntos:** si las filas tienen soportes que no se
    superponen, el producto interno es cero porque no existe ninguna
    posición en la que ambas tengan valores no nulos.
2.  **Soportes anidados:** si el soporte de una fila está contenido en
    una mitad constante de otra fila de mayor escala, el producto
    interno también es cero. La fila de mayor escala tiene el mismo
    signo en el soporte menor, mientras que la fila de menor escala
    contiene la misma cantidad de valores positivos y negativos; por
    ello, la suma se cancela.

También la fila de aproximación global es ortogonal a cualquier fila de
detalle, porque la suma de los componentes positivos y negativos de cada
detalle es cero.

Por ejemplo, para las dos primeras filas de $H_8$, el producto interno
es:

```math
\begin{aligned} \langle h_0,h_1\rangle &= 4\left(\frac{1}{\sqrt{8}}\frac{1}{\sqrt{8}}\right) + 4\left(\frac{1}{\sqrt{8}}\frac{-1}{\sqrt{8}}\right)\\ &=\frac{4}{8}-\frac{4}{8}=0. \end{aligned}
```

Por tanto, las filas distintas son ortogonales. Al combinar la
ortogonalidad con la norma unitaria se obtiene:

```math
\boxed{H_NH_N^T=I_N}.
```

Esta propiedad permite utilizar la transpuesta como matriz inversa y
garantiza que la transformación sea reversible cuando se conservan todos
los coeficientes.

## 4. Transformada de Haar bidimensional

Una imagen digital en escala de grises se representa mediante una
matriz:

```math
A\in\mathbb{R}^{M\times N},
```

donde cada elemento $a_{ij}$ representa la intensidad del píxel ubicado
en la fila $i$ y columna $j$.

La extensión de Haar a dos dimensiones se obtiene aplicando la
transformación unidimensional sobre ambos ejes. Esta propiedad se
denomina separabilidad: la transformación bidimensional puede calcularse
mediante dos productos matriciales consecutivos \[3\], \[4\].

Sean $H_M$ la matriz de Haar de orden $M$ y $H_N$ la matriz de Haar de
orden $N$. La transformación directa se define como:

```math
C=H_MAH_N^T.
```

La operación puede separarse en dos etapas:

```math
B=H_MA,
```

```math
C=BH_N^T.
```

En la primera etapa se transforman las filas de la matriz de imagen,
según la convención de multiplicación adoptada. En la segunda etapa se
transforman las columnas. El resultado $C$ contiene los coeficientes de
aproximación y detalle en los dos ejes.

Si la matriz de imagen es cuadrada y ambas dimensiones tienen el mismo
tamaño, se utiliza la misma matriz de transformación en ambos ejes:

```math
C=H_NAH_N^T.
```

### 4.1. Transformación inversa bidimensional

Si las matrices $H_M$ y $H_N$ son ortonormales, se cumple:

```math
H_M^{-1}=H_M^T, \qquad H_N^{-1}=H_N^T.
```

Partiendo de la transformación directa:

```math
C=H_MAH_N^T,
```

se multiplica por $H_M^T$ a la izquierda y por $H_N$ a la derecha:

```math
H_M^TCH_N = H_M^TH_MA H_N^TH_N.
```

Como:

```math
H_M^TH_M=I_M, \qquad H_N^TH_N=I_N,
```

se obtiene:

```math
\boxed{A=H_M^TCH_N}.
```

Esta expresión demuestra que la transformación bidimensional es
reversible si se mantienen todos los coeficientes.

## 5. Interpretación matemática de las subbandas

En una descomposición bidimensional de un nivel, cada eje se divide en
componentes de aproximación $L$ y de detalle $H$. Al combinar las dos
posibilidades de cada eje aparecen cuatro subbandas:

-   **LL:** aproximación en ambos ejes. Representa la componente de baja
    variación espacial.
-   **LH:** aproximación en un eje y detalle en el otro.
-   **HL:** detalle en un eje y aproximación en el otro.
-   **HH:** detalle en ambos ejes.

Estas subbandas corresponden a productos tensoriales de las funciones de
escala y wavelet. Si $\phi$ representa la función de escala y
$\psi$ la wavelet madre, las cuatro familias bidimensionales
pueden expresarse como:

```math
\phi(x)\phi(y),\qquad \phi(x)\psi(y),\qquad \psi(x)\phi(y),\qquad \psi(x)\psi(y).
```

La primera combinación representa aproximación en ambas coordenadas; las
otras tres representan variaciones o detalles en uno o en ambos ejes. La
denominación horizontal o vertical de LH y HL depende de qué eje se
transforma primero y de la convención de coordenadas utilizada. Por
ello, al describir una implementación se debe declarar explícitamente
esa convención \[4\].

En una descomposición multinivel, el procedimiento se aplica nuevamente
sobre la subbanda LL. Así, la aproximación obtenida en un nivel se
divide en una nueva aproximación y tres nuevos grupos de detalles. Este
proceso genera una representación jerárquica de la información \[1\],
\[3\].

## 6. Análisis multirresolución

El análisis multirresolución proporciona el marco matemático para
interpretar la descomposición wavelet. Consiste en representar una
función mediante espacios de aproximación anidados:

```math
\cdots\subset V_{-1}\subset V_0\subset V_1\subset\cdots
```

Cada espacio $V_j$ contiene aproximaciones de una función a una
determinada resolución. El espacio $W_j$ representa la información
adicional necesaria para pasar de una resolución a la siguiente. En el
análisis ortogonal se cumple:

```math
V_{j+1}=V_j\oplus W_j.
```

El símbolo $\oplus$ indica suma directa ortogonal: los
componentes de $V_j$ y $W_j$ son ortogonales entre sí y juntos generan
$V_{j+1}$.

En la construcción de Haar, la función de escala genera los espacios de
aproximación y la wavelet madre genera los espacios de detalle. Las
versiones trasladadas y dilatadas de estas funciones se definen como:

```math
\phi_{j,k}(x)=2^{j/2}\phi(2^jx-k),
```

```math
\psi_{j,k}(x)=2^{j/2}\psi(2^jx-k).
```

Aquí, $j$ identifica la escala y $k$ identifica la traslación. El factor
$2^{j/2}$ normaliza las funciones para conservar la norma unitaria.

La descomposición permite expresar una función mediante una aproximación
en una escala determinada y los detalles correspondientes a las escalas
más finas. Este principio constituye la base teórica de la
representación multinivel de Haar \[1\].

## 7. Conservación de energía e identidad de Parseval

Una transformación ortonormal conserva la norma euclidiana de los
vectores. Si $x$ es el vector original y $c=H_Nx$ su representación en
la base de Haar, entonces:

```math
\|c\|_2^2 = c^Tc = (H_Nx)^T(H_Nx).
```

Desarrollando:

```math
\|c\|_2^2 = x^TH_N^TH_Nx.
```

Como $H_N^TH_N=I_N$, resulta:

```math
\boxed{\|c\|_2^2=\|x\|_2^2}.
```

Esta igualdad corresponde a la identidad de Parseval: la suma de los
cuadrados de los valores originales es igual a la suma de los cuadrados
de sus coeficientes transformados.

Para matrices, la misma propiedad se expresa mediante la norma de
Frobenius. Si $C=H_MAH_N^T$, entonces:

```math
\boxed{\|C\|_F^2=\|A\|_F^2}.
```

Por tanto, la transformación redistribuye la energía entre coeficientes,
pero no la modifica mientras se mantengan todos ellos. Esta propiedad es
consecuencia directa de la ortonormalidad \[1\], \[3\].

## 8. Propiedad de reconstrucción y error en el dominio transformado

La transformación completa es reversible. Sin embargo, si durante un
proceso posterior algunos coeficientes se modifican o se sustituyen por
cero, la reconstrucción ya no tiene por qué coincidir con la matriz
original.

Sea $C$ la matriz completa de coeficientes y $\widetilde{C}$ una
matriz modificada. La reconstrucción correspondiente es:

```math
\widetilde A=H_M^T\widetilde C H_N.
```

La diferencia entre la imagen original y la reconstruida es:

```math
A-\widetilde A = H_M^T(C-\widetilde C)H_N.
```

Como las matrices de transformación son ortonormales, la norma de
Frobenius se conserva. Por tanto:

```math
\boxed{ \|A-\widetilde A\|_F^2 = \|C-\widetilde C\|_F^2 }.
```

Esta igualdad establece que el error cuadrático de reconstrucción en el
dominio original es exactamente igual al error cuadrático entre los
coeficientes originales y los modificados en el dominio transformado. Si
los coeficientes descartados se reemplazan por cero, el error
corresponde a la suma de sus cuadrados.

La propiedad explica matemáticamente por qué la magnitud de los
coeficientes es relevante al estudiar una selección basada en error
cuadrático. No obstante, la transformación de Haar por sí sola no
determina un formato de archivo comprimido: la reducción efectiva del
almacenamiento también depende de la cuantización, la codificación y la
información necesaria para identificar los coeficientes conservados.

## 9. Referencias bibliográficas

\[1\] S. G. Mallat, "Multiresolution Approximations and Wavelet Orthonormal Bases of $L^2(\mathbb{R})$," *Transactions of the American Mathematical Society*, vol. 315, no. 1, pp. 69--87, 1989, doi: 10.1090/S0002-9947-1989-1008470-5.

\[2\] A. Haar, "Zur Theorie der orthogonalen Funktionensysteme," *Mathematische Annalen*, vol. 69, pp. 331--371, 1910, doi: 10.1007/BF01456326.

\[3\] I. Daubechies, *Ten Lectures on Wavelets*. Philadelphia, PA, USA: Society for Industrial and Applied Mathematics, 1992, doi: 10.1137/1.9781611970104.

\[4\] PyWavelets Developers, "2D Forward and Inverse Discrete Wavelet Transform," *PyWavelets Documentation*. \[En línea\]. Disponible en: https://pywavelets.readthedocs.io/en/stable/ref/2d-dwt-and-idwt.html. \[Consultado: 17-nov-2026\].
