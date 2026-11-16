# Transformada Discreta del Coseno Tipo II (DCT-II)



## 1\. Introducción



En el análisis de señales y el procesamiento digital de imágenes, las transformadas lineales cumplen el propósito fundamental de reducir la complejidad de un problema al proyectar los datos hacia un dominio alternativo donde sus

componentes fundamentales resultan más fáciles de manipular. Si bien es frecuente asumir que la Transformada Discreta del Coseno (DCT) es simplemente una aproximación discreta de la Transformada Coseno de Fourier continua (FCT), dicha interpretación es matemáticamente imprecisa. Como demostró Gilbert Strang, la familia completa de transformadas trigonométricas discretas (DCT y DST) emerge de manera natural como las soluciones propias de un sistema físico discretizado: el oscilador armónico simple sometido a distintas combinaciones de condiciones de frontera homogéneas.   El presente informe desarrolla formalmente el origen diferencial y matricial de la DCT Tipo II (DCT-II), establece sus diferencias estructurales frente a otras transformadas ortogonales, demuestra su propiedad de ortonormalidad y expone las factorizaciones de matrices dispersas que permiten su cálculo rápido en bloques bidimensionales de $8 \\times 8$ para la compresión digital



## 2\. Origen Fisico-Mátematico

#### **2.1 Del oscilador Armónico a la Matriz $A\_2$**

Para deducir la DCT-II desde primeros principios, se parte del problema continuo de autovalores asociado al oscilador armónico simple no amortiguado en el dominio espacial $x \\in \[0, \\pi]$:



&#x09;				$$\\frac{d^2u}{dx^2} + \\lambda u = 0 \\quad \\Longleftrightarrow \\quad u'' + \\lambda u = 0$$



Al imponer condiciones de frontera homogéneas de Neumann (derivada nula en ambos extremos, $u'(0) = 0$ y $u'(\\pi) = 0$), las autofunciones del sistema continuo resultan ser los modos naturales de vibración $\\cos(kx)$ con autovalores $\\lambda = k^2$ para $k = 0, 1, 2, \\dots$.   Al trasladar este modelo al dominio discreto sobre una malla equiespaciada de $N$ puntos indexados por $l = 0, 1, \\dots, N-1$, la segunda derivada se aproxima mediante diferencias finitas centrales:

&#x09;			$$-u\_{l-1} + 2u\_l - u\_{l+1} = \\lambda u\_l$$



La distinción entre las ocho variantes de la DCT radica en el punto exacto donde se imponen las condiciones de frontera de simetría (Neumann) o antisimetría (Dirichlet). Específicamente, la DCT-II surge cuando las condiciones de frontera de Neumann se aplican en los puntos medios de la malla (mid-grid) en ambos extremos, es decir, en $l = -1/2$ y $l = N - 1/2$. Esto equivale a imponer una simetría de medio paso (Half-Sample Symmetric o HSHS) dada por $u\_{-1} = u\_0$ en el extremo izquierdo y $u\_N = u\_{N-1}$ en el extremo derecho.



Al sustituir estas condiciones de borde en la ecuación en diferencias para $l = 0$ y $l = N-1$, el sistema de ecuaciones lineales se expresa en forma matricial como un problema de autovalores $A\_2 \\vec{u}\_k = \\lambda\_k \\vec{u}\_k$, donde $A\_2$ es una matriz tridiagonal simétrica de diferencias de segundo orden de orden $N \\times N$:



$$A\_2 = \\begin{pmatrix} 1 \& -1 \& 0 \& \\cdots \& 0 \\\\ -1 \& 2 \& -1 \& \\cdots \& 0 \\\\ 0 \& -1 \& 2 \& \\ddots \& \\vdots \\\\ \\vdots \& \\vdots \& \\ddots \& 2 \& -1 \\\\ 0 \& 0 \& \\cdots \& -1 \& 1 \\end{pmatrix}$$



Resolviendo este sistema algebraico, se verifica que los autovalores de $A\_2$ están dados por $\\lambda\_k = 2 - 2\\cos\\left(\\frac{k\\pi}{N}\\right)$ y que la $l$-ésima componente del $k$-ésimo autovector $\\vec{u}\_k$ corresponde exactamente al núcleo cosenoidal de la DCT-II:  



&#x09;	 $$(u\_k)\_l = \\cos\\left\[\\left(l + \\frac{1}{2}\\right)\\frac{k\\pi}{N}\\right], \\quad l, k = 0, 1, \\dots, N-1$$





## 3\. Diferencias Estructurales con Otras Transformadas Ortogonales



La elección de la DCT-II frente a otras transformadas lineales se fundamenta en cuatro propiedades distintivas:



1. **Frente a la Transformada Discreta de Fourier (DFT):** La DFT asume que una señal finita se repite de manera periódica circular. Si los valores inicial y final de un bloque de datos difieren, la extensión circular introduce un salto discontinuo artificial en los bordes, dispersando energía hacia frecuencias altas (fenómeno de Gibbs). Por el contrario, la extensión simétrica en punto medio (HSHS) de la DCT-II garantiza continuidad en los bordes del bloque, concentrando la energía en unos pocos coeficientes de baja frecuencia. 

2. **Frente a las demás variantes de la familia DCT/DST:** Mientras la DCT-I impone simetría de Neumann directamente sobre las muestras extremas (WSWS), las variantes DCT-III y DCT-IV reemplazan la condición del extremo derecho por una condición de Dirichlet ($u(\\pi)=0$, antisimetría). Por su parte, la familia de Transformadas Discretas del Seno (DST-I a DST-VIII) sustituye la condición de Neumann en el origen por una condición de Dirichlet. 

3. **Frente a transformadas de funciones cuadradas (Walsh-Hadamard y Haar):** A diferencia de la Transformada de Walsh-Hadamard (WHT) o la de Haar —cuyas funciones base son pulsos escalonados binarios compuestos únicamente por valores $\\{+1, -1, 0\\}$ provenientes del álgebra combinatoria y el análisis multirresolución—, la DCT-II utiliza funciones base sinusoidales suaves. Esto le permite aproximar señales naturales continuas con un error cuadrático medio significativamente menor.  

4. **Aproximación óptima a la Transformada de Karhunen-Loéve (KLT):** Aunque la KLT es teóricamente óptima para descorrelacionar cualquier señal estadística, su matriz de transformación depende de la covarianza de cada imagen y carece de algoritmos rápidos. La DCT-II utiliza una base fija independiente de los datos que converge asintóticamente a la KLT para procesos estocásticos de Markov de primer orden (modelo estándar de las imágenes naturales), logrando una compactación de energía casi idéntica.   



## 4\. Definición matricial y ortonormalidad

Para garantizar que la transformación conserve la energía de la señal (**Teorema de Parseval**) y sea fácilmente invertible, los autovectores $\\vec{u}\_k$ de la matriz $A\_2$ deben normalizarse. De este modo, los elementos $(k, n)$ de la matriz ortonormal de la DCT-II de orden $N$, denotada como $C\_N^{II}$, se definen como:



&#x20;            $$\[C\_N^{II}]\_{kn} = \\sqrt{\\frac{2}{N}} \\epsilon\_k \\cos\\left\[\\frac{\\pi(2n+1)k}{2N}\\right], \\quad k, n = 0, 1, \\dots, N-1$$



donde el factor de normalización $\\epsilon\_k$ está definido por:



&#x09;	$$\\epsilon\_k = \\begin{cases} \\dfrac{1}{\\sqrt{2}} \& \\text{si } k = 0 \\\\ 1 \& \\text{si } k \\in \\{1, 2, \\dots, N-1\\} \\end{cases}$$



Para demostrar que $C\_N^{II}$ es una matriz ortogonal (unitaria real), debe probarse que el producto interno entre cualesquiera dos vectores fila $\\vec{\\psi}\_k$ y $\\vec{\\psi}\_m$ satisface $\\langle \\vec{\\psi}\_k, \\vec{\\psi}\_m \\rangle = \\delta\_{km}$.



**1. Ortogonalidad ($k \\neq m$):**

De acuerdo con la deducción de la Sección 2, las filas no escaladas de $C\_N^{II}$ son los autovectores de la matriz $A\_2$. Dado que $A\_2$ es una matriz real y simétrica ($A\_2 = A\_2^T$), el Teorema Espectral del álgebra lineal garantiza que los autovectores correspondientes a autovalores distintos ($\\lambda\_k \\neq \\lambda\_m$ para $k \\neq m$) son mutuamente ortogonales:



&#x09;	$$\\sum\_{n=0}^{N-1} \\cos\\left\[\\frac{\\pi(2n+1)k}{2N}\\right] \\cos\\left\[\\frac{\\pi(2n+1)m}{2N}\\right] = 0 \\quad \\forall \\, k \\neq m$$



**2. Normalización ($k = m$):**

Evaluando la norma al cuadrado de cada fila:

Para $k = m = 0$: Como $\\cos(0) = 1$ y $\\epsilon\_0 = \\frac{1}{\\sqrt{2}}$, se obtiene:

&#x09;	

&#x09;	$$\\langle \\vec{\\psi}\_0, \\vec{\\psi}\_0 \\rangle = \\sum\_{n=0}^{N-1} \\left(\\sqrt{\\frac{2}{N}} \\cdot \\frac{1}{\\sqrt{2}} \\cdot 1\\right)^2 = \\sum\_{n=0}^{N-1} \\frac{1}{N} = 1$$



Para $k = m > 0$: Aplicando la identidad trigonométrica $\\cos^2(\\theta) = \\frac{1 + \\cos(2\\theta)}{2}$ y $\\epsilon\_k = 1$:

&#x09;

&#x09;	$$\\langle \\vec{\\psi}\_k, \\vec{\\psi}\_k \\rangle = \\frac{2}{N} \\sum\_{n=0}^{N-1} \\cos^2\\left\[\\frac{\\pi(2n+1)k}{2N}\\right] = \\frac{1}{N} \\sum\_{n=0}^{N-1} \\left(1 + \\cos\\left\[\\frac{\\pi(2n+1)k}{N}\\right]\\right)$$



Dado que la suma de una onda cosenoidal completa sobre un número entero de semiperiodos se anula:



&#x09;	$\\sum\_{n=0}^{N-1} \\cos\\left\[\\frac{\\pi(2n+1)k}{N}\\right] = 0$ 



para $0 < k < N$), la expresión se reduce a:



&#x09;		$$\\langle \\vec{\\psi}\_k, \\vec{\\psi}\_k \\rangle = \\frac{1}{N} (N + 0) = 1$$



Queda demostrado que $C\_N^{II} (C\_N^{II})^T = I\_N$. En consecuencia, la matriz inversa de la DCT-II coincide exactamente con su matriz transpuesta, la cual define a la DCT Tipo III (DCT-III): 







## 5\. Algoritmos Rápidos de Cálculo (Fast DCT-II)



El cálculo directo de la transformada unidimensional $\\vec{c} = C\_N^{II} \\vec{x}$ mediante multiplicación matriz-vector convencional requiere $N^2$ multiplicaciones y $N(N-1)$ sumas. Para reducir esta complejidad computacional de $\\mathcal{O}(N^2)$ a $\\mathcal{O}(N \\log\_2 N)$ cuando $N = 2^m$, los algoritmos rápidos aprovechan la simetría par e impar de los vectores base:  



&#x09;	$$\\Psi\_{k, N-1-n} = (-1)^k \\Psi\_{k, n}, \\quad n = 0, 1, \\dots, \\frac{N}{2}-1$$



Esta propiedad indica que las filas de índice par ($2k$) son simétricas respecto al centro del bloque, mientras que las filas de índice impar ($2k+1$) son antisimétricas.



#### **5.1 Factorización Recursiva en Matrices Dispersas**

Separando la sumatoria de la DCT-II (omitiendo temporalmente los factores de normalización) en coeficientes pares e impares, y agrupando los términos extremos de la señal de entrada $x\_n$ y $x\_{N-1-n}$, el cálculo de tamaño $N$ se divide en dos subproblemas independientes de tamaño $N/2$:



**Coeficientes de índice par ($2k$):**



&#x09;	$$c\_{2k}^{II} = \\sum\_{n=0}^{\\frac{N}{2}-1} (x\_n + x\_{N-1-n}) \\cos\\left\[\\frac{\\pi(2n+1)k}{2(N/2)}\\right], \\quad k = 0, 1, \\dots, \\frac{N}{2}-1$$



Esta expresión corresponde exactamente a una DCT-II de longitud $N/2$ aplicada sobre la suma de las muestras simétricas.



**Coeficientes de índice impar ($2k+1$):**



$$c\_{2k+1}^{II} = \\sum\_{n=0}^{\\frac{N}{2}-1} (x\_n - x\_{N-1-n}) \\cos\\left\[\\frac{\\pi(2n+1)(2k+1)}{4(N/2)}\\right], \\quad k = 0, 1, \\dots, \\frac{N}{2}-1$$



Esta expresión corresponde a una DCT Tipo IV (DCT-IV) de longitud $N/2$ aplicada sobre la diferencia de las muestras simétricas.



En notación matricial, esta descomposición recursiva se expresa mediante el producto de matrices dispersas (sparse matrices).





$$C\_N^{II} = P\_N \\begin{pmatrix} C\_{N/2}^{II} \& 0 \\\\ 0 \& J\_{N/2} C\_{N/2}^{IV} J\_{N/2} \\end{pmatrix} \\begin{pmatrix} I\_{N/2} \& J\_{N/2} \\\\ J\_{N/2} \& -I\_{N/2} \\end{pmatrix}$$





donde $I\_{N/2}$ es la matriz identidad, $J\_{N/2}$ es la matriz de reflexión antidiagonal (que invierte el orden de los elementos) y $P\_N$ es una matriz de permutación que reordena los coeficientes pares e impares.



#### **5.2 Otras Estrategias Algorítmicas y Cota Óptima para $N=8$**



Además de la factorización directa anterior (algoritmo de Chen), el desarrollo algorítmico permite implementar la DCT-II mediante tres enfoques alternativos de alta eficiencia:



1. **Algoritmos Split-Radix:** Descomponen la DCT-II de $N$ puntos en una DCT-II de $N/2$, una DCT-II de $N/4$ y una DST-II de $N/4$, minimizando las operaciones redundantes.

2. **Cálculo vía Transformada de Walsh-Hadamard (WHT):** Dado que la WHT posee la misma estructura de simetría par/impar que la DCT-II, se aplica primero la WHT (que solo utiliza sumas y restas con elementos $\\pm 1$) y posteriormente una matriz de conversión dispersa por bloques con rotaciones de plano**.**

3. **Cálculo vía DFT Real (FFCT):** Reordena la secuencia de entrada para mapear la DCT-II sobre una Transformada Rápida de Fourier (FFT) real seguida de rotaciones angulares.



Para el caso estándar de $N = 8$, la cota teórica mínima de multiplicaciones no triviales está determinada por $\\mu(C\_N^{II}) = 2^{n+1} - n - 2$ con $N = 2^n$. Al factorizar los factores de escala $\\sqrt{2}$ hacia el final del proceso (DCT-II escalada, como en el algoritmo de Loeffler), se alcanza el mínimo teórico absoluto de 11 multiplicaciones y 29 sumas.



## 6\. Extensión Bidimensional ($8 \\times 8$) y Aplicación a la Compresión Digital



En el procesamiento de imágenes y codificación de video (estándares como JPEG, MPEG y H.264), los datos se procesan particionando la imagen en bloques bidimensionales de tamaño $N \\times N$, siendo $8 \\times 8$ píxeles el estándar universal.



La elección de bloques de $8 \\times 8$ ($64$ píxeles) responde a un compromiso óptimo entre estadística de la señal y complejidad computacional:



* **Potencia de dos ($N = 2^3$):** Permite aplicar tres niveles completos de descomposición recursiva de matrices dispersas hasta llegar a operaciones triviales de $2$ puntos.   



* **Correlación espacial vs. Artefactos:** En bloques menores ($4 \\times 4$), la compactación de energía es insuficiente; en bloques mayores ($16 \\times 16$ o superiores), la presencia de bordes abruptos rompe la correlación local de Markov-1 entre píxeles distantes, introduciendo ruido alrededor de los contornos (ringing) e incrementando el costo aritmético.



#### **6.1 Formulación Matricial y Producto de Kronecker**



Gracias a la separabilidad del núcleo bidimensional de la DCT-II, la transformación de una matriz de píxeles $X \\in \\mathbb{R}^{8 \\times 8}$ en una matriz de coeficientes espectrales $Y \\in \\mathbb{R}^{8 \\times 8}$ se realiza aplicando la transformada 1-D primero sobre las filas y luego sobre las columnas: 



&#x09;	$$Y = C\_8^{II} X (C\_8^{II})^T$$



Si se vectorizan las matrices $X$ e $Y$ por concatenación lexicográfica de sus filas, el operador bidimensional se expresa de forma compacta mediante el producto de Kronecker ($\\otimes$):



&#x09;	$$\\text{vec}(Y) = (C\_8^{II} \\otimes C\_8^{II}) \\text{vec}(X)$$



Mientras el método directo requiere $64 \\times 64 = 4096$ multiplicaciones por bloque y el método convencional fila-columna requiere $16 \\times 11 = 176$ multiplicaciones, el enfoque estructural basado en la factorización del producto de Kronecker y rotaciones de Givens-Jacobi alcanza la cota teórica mínima de 88 multiplicaciones para todo el bloque $8 \\times 8$.



#### **6.2 Proceso de Comprensión**



La transformación de una imagen es una clase de matrices unitarias utilizadas para representar imágenes. Así como una señal 1D se expresa mediante funciones base ortogonales, una imagen puede expandirse usando un conjunto discreto de arreglos llamados imágenes base *(basis images)*.



**Vectorialmente:**

Una imagen de tamaño $N \\times N$ también puede verse como un vector de dimensión $N^2 \\times 1$, donde la transformada proporciona un nuevo sistema de coordenadas (vectores base) para ese espacio vectorial.



**Transformación unitaria: 1D**

Para una secuencia unidimensional $u(n)$ representada como un vector $\\mathbf{u}$ de tamaño $N$:

* Transformada directa:



&#x09;	$\\mathbf{v} = \\mathbf{A}\\mathbf{u} \\implies v(k) = \\sum\_{n=0}^{N-1} a(k,n)u(n)$. 



* Transformada inversa: Dado que la matriz $\\mathbf{A}$ es unitaria, su inversa es igual a su transpuesta conjugada ($\\mathbf{A}^{-1} = \\mathbf{A}^{\*T}$). Por lo tanto, la señal original se recupera como:



&#x09;	$$\\mathbf{u} = \\mathbf{A}^{\*T}\\mathbf{v} \\implies u(n) = \\sum\_{k=0}^{N-1} v(k)a^\*(k,n)$$



**Vectores base:**

Las columnas de la matriz inversa $\\mathbf{A}^{\*T}$ (denotadas como $\\mathbf{a}\_k^\*$) se llaman vectores base (basis vectors) de $\\mathbf{A}$. La ecuación inversa (5.2) muestra que cualquier secuencia $u(n)$ se puede reconstruir como una suma ponderada de estos vectores base. 

Utilidad de los coeficientes $v(k)$: Los coeficientes obtenidos representan la señal original en otro dominio y son fundamentales para tareas como filtrado, compresión de datos y extracción de características.



**Tranformación Bidimensaional 2D**

Para una imagen bidimensional $u(m,n)$ de tamaño $N \\times N$, la expansión ortogonal general se define mediante dos ecuaciones: 

* **Transformada Directa:**

Convierte los píxeles $u(m,n)$ en coeficientes de transformada $v(k,l)$, los cuales forman la matriz $\\mathbf{V}$ llamada imagen transformada:



&#x09;	$$v(k,l) = \\sum\_{m,n=0}^{N-1} u(m,n) \\, a\_{k,l}(m,n), \\quad 0 \\le k, l \\le N-1$$



* **Transformada Inversa:**

Reconstruye la imagen original a partir de los coeficientes:



&#x09;	$$u(m,n) = \\sum\_{k,l=0}^{N-1} v(k,l) \\, a^\*\_{k,l}(m,n), \\quad 0 \\le m, n \\le N-1$$



**Propiedades del Núcleo de Transformación $\\{a\_{k,l}(m,n)\\}$:**



El conjunto de funciones base discretas debe cumplir dos condiciones fundamentales:



1. **Ortonormalidad:** Las funciones base son perpendiculares entre sí y tienen energía unitaria.



&#x09;		$\\sum a\_{k,l}(m,n) a^\*\_{k',l'}(m,n) = \\delta(k-k', l-l')$



Esto evita redundancia entre los coeficientes.



**2.   Completitud:** El conjunto de funciones base es suficiente para representar cualquier imagen de tamaño $N \\times N$ sin que falte información.



&#x09;	$\\sum a\_{k,l}(m,n) a^\*\_{k,l}(m',n') = \\delta(m-m', n-n')$



Importancia para la Compresión (Truncamiento y Error Mínimo):



**Expansión Truncada:** : Si en lugar de usar todos los $N \\times N$ coeficientes solo usamos un subconjunto más pequeño ($P \\le N$ y $Q \\le N$) para aproximar la imagen ($u\_{P,Q}(m,n)$), la propiedad de ortonormalidad garantiza que los coeficientes $v(k,l)$ calculados en (5.3) minimizan la suma del error cuadrático.



&#x09;	$$\\sigma\_e^2 = \\sum\_{m,n=0}^{N-1} \[u(m,n) - u\_{P,Q}(m,n)]^2$$



**Reconstrucción Exacta:** La propiedad de completitud asegura que si no descartamos ningún coeficiente ($P = Q = N$), el error $\\sigma\_e^2$ será exactamente cero (compresión sin pérdida de información en la transformación).



**Separabilidad:**



Calcular directamente una transformada 2D general (ecuación 5.3) requiere un número de multiplicaciones y sumas del orden de $O(N^4)$, lo cual resulta excesivo e impráctico para imágenes de tamaño real.  La complejidad computacional se reduce drásticamente de $O(N^4)$ a $O(N^3)$ (e incluso menos si se usan algoritmos rápidos) cuando se restringe la transformada para que sea separable.

Un núcleo 2D $a\_{k,l}(m,n)$ es separable si puede expresarse como el producto de dos funciones base unidimensionales (una para las filas y otra para las columnas):



&#x09;	$$a\_{k,l}(m,n) = a\_k(m) \\, b\_l(n) \\triangleq a(k,m) \\, b(l,n)$$



**Matrices Unitarias 1D:**



Para que la transformada 2D siga siendo ortonormal y completa (cumpliendo las ecuaciones 5.5 y 5.6), las matrices unidimensionales $\\mathbf{A} \\triangleq \\{a(k,m)\\}$ y $\\mathbf{B} \\triangleq \\{b(l,n)\\}$ deben ser matrices unitarias por sí mismas, es decir:



&#x09;	$$\\mathbf{A}\\mathbf{A}^{\*T} = \\mathbf{A}^T\\mathbf{A}^\* = \\mathbf{I}$$



Si $\\mathbf{A}$ y $\\mathbf{B}$ son unitarias individualmente, sus vectores fila/columna ya son ortonormales y completos en 1D, lo cual garantiza automáticamente que el producto en 2D también lo sea.



Por lo general, se elige la misma matriz de transformación para ambas dimensiones ($\\mathbf{B} = \\mathbf{A}$), lo que permite calcular la transformada 2D mediante una simple multiplicación de matrices ($\\mathbf{V} = \\mathbf{A}\\mathbf{U}\\mathbf{A}^T$).



#### 6.1 Imagen cuadrada ($N \\times N$):

* **Transformada directa:**



&#x09;	$$v(k,l) = \\sum\_{m,n=0}^{N-1} a(k,m) u(m,n) a(l,n) \\iff \\mathbf{V} = \\mathbf{A}\\mathbf{U}\\mathbf{A}^T$$



* **Transformada inversa:**


&#x09;	$$u(m,n) = \\sum\_{k,l=0}^{N-1} a^\*(k,m) v(k,l) a^\*(l,n) \\iff \\mathbf{U} = \\mathbf{A}^{\*T}\\mathbf{V}\\mathbf{A}^\*$$



#### 6.2 Imagen rectangular ($M \\times N$):

Se usan dos matrices unitarias de distinto tamaño ($\\mathbf{A}\_M$ de $M\\times M$ y $\\mathbf{A}\_N$ de $N\\times N$): $\\mathbf{V} = \\mathbf{A}\_M \\mathbf{U} \\mathbf{A}\_N$ (5.13) y su inversa $\\mathbf{U} = \\mathbf{A}\_M^{\*T} \\mathbf{V} \\mathbf{A}\_N^{\*T}$.

* **Transformada directa:**



&#x09;	$$\\mathbf{V} = \\mathbf{A}\_M \\mathbf{U} \\mathbf{A}\_N$$



* **Transformada inversa**



&#x09;	$$\\mathbf{U} = \\mathbf{A}\_M^{\*T} \\mathbf{V} \\mathbf{A}\_N^{\*T}$$



**Calculo:**

Al reescribir la directa como $\\mathbf{V}^T = \\mathbf{A}\[\\mathbf{A}\\mathbf{U}]^T$, se demuestra que una transformada 2D se ejecuta aplicando primero la transformada 1D a cada columna de $\\mathbf{U}$ ($\\mathbf{A}\\mathbf{U}$) y luego aplicando la transformada 1D a cada fila del resultado para obtener $\\mathbf{V}$.



&#x09;	$$\\mathbf{V}^T = \\mathbf{A}\[\\mathbf{A}\\mathbf{U}]^T$$





Para una transformación bidimensional separable, las matrices base fundamentales se construyen a partir de los vectores base unidimensionales $\\mathbf{a}\_k^\*$ (correspondientes a las columnas de $\\mathbf{A}^{\*T}$). Mediante el producto exterior de estos vectores se define el conjunto completo de $N^2$ imágenes base ($\\mathbf{A}\_{k,l}^\*$) de dimensión $N \\times N$:  



&#x09;	$$\\mathbf{A}\_{k,l}^\* = \\mathbf{a}\_k^\* \\mathbf{a}\_l^{\*T}$$



Con el fin de establecer una representación en serie para una imagen $\\mathbf{U}$, se define el producto interno matricial entre dos arreglos $\\mathbf{F}$ y $\\mathbf{G}$ de tamaño $N \\times N$ como:



&#x09;	$$\\langle \\mathbf{F}, \\mathbf{G} \\rangle = \\sum\_{m=0}^{N-1} \\sum\_{n=0}^{N-1} f(m,n) g^\*(m,n)$$



A partir de este operador, cualquier imagen $\\mathbf{U}$ puede expresarse como una combinación lineal ponderada de las $N^2$ imágenes base ortonormales:   



&#x09;	$$\\mathbf{U} = \\sum\_{k,l=0}^{N-1} v(k,l) \\mathbf{A}\_{k,l}^\*$$



En esta expansión, cada coeficiente de la transformada $v(k,l)$ representa geométricamente la proyección de la imagen original $\\mathbf{U}$ sobre la $(k,l)$-ésima imagen base, calculada mediante su producto interno:



&#x09;	$$v(k,l) = \\langle \\mathbf{U}, \\mathbf{A}\_{k,l}^\* \\rangle$$



Adicionalmente, si las matrices $\\mathbf{U}$ y $\\mathbf{V}$ se reordenan por filas en vectores unidimensionales $\\mathbf{u}$ y $\\mathbf{v}$ de dimensión $N^2 \\times 1$, el par de transformadas bidimensionales adopta una formulación lineal compacta utilizando el producto de Kronecker ($\\otimes$):



&#x09;	$$\\mathbf{v} = (\\mathbf{A} \\otimes \\mathbf{A})\\mathbf{u} \\triangleq \\mathcal{A}\\mathbf{u}$$



&#x09;	$$\\mathbf{u} = (\\mathbf{A} \\otimes \\mathbf{A})^{\*T}\\mathbf{v} = \\mathcal{A}^{\*T}\\mathbf{v}$$



Donde $\\mathcal{A} \\triangleq \\mathbf{A} \\otimes \\mathbf{A}$ representa el operador unitario global en el espacio vectorial de dimensión $N^2$.



#### 6.3 Separabilidad mediante Producto de Kronecker

Una transformación lineal arbitraria $\\boldsymbol{y} = \\mathcal{A}\\boldsymbol{x}$ (5.23) se denomina separable si su matriz $\\mathcal{A}$ (de tamaño $N^2 \\times N^2$) puede descomponerse como el producto de Kronecker de dos matrices más pequeñas ($\\mathbf{A}\_1$ y $\\mathbf{A}\_2$ de tamaño $N \\times N$):



&#x09;	$$\\mathcal{A} = \\mathbf{A}\_1 \\otimes \\mathbf{A}\_2$$



Esto permite reducir el cálculo vectorial a una transformación bidimensional matricial mucho más eficiente:



&#x09;	$$\\mathbf{Y} = \\mathbf{A}\_1 \\mathbf{X} \\mathbf{A}\_2^T$$



Donde $\\mathbf{X}$ y $\\mathbf{Y}$ son las matrices que al ordenarse por filas forman los vectores $\\boldsymbol{x}$ y $\\boldsymbol{y}$).



**Dimensionalidad y Transformadas Rápidas**

Los $2N^3$ cálculos necesarios para obtener $\\mathbf{V}$ se pueden reducir aún más utilizando transformadas rápidas (fast transforms), cuya matriz $\\mathbf{A}$ se factoriza en un producto de $p$ matrices dispersas (con muy pocos elementos distintos de cero, $r \\ll N$):



&#x09;	$$\\mathbf{A} = \\mathbf{A}\_{(1)} \\mathbf{A}\_{(2)} \\dots \\mathbf{A}\_{(p)}$$



Esto reduce las operaciones de una multiplicación $\\mathbf{y} = \\mathbf{A}\\mathbf{x}$ a $rpN$. Para transformadas como la DCT (coseno), Fourier, seno, Hadamard y Slant, donde $p = \\log\_2 N$, la complejidad baja al orden de $N \\log\_2 N$ (o $N^2 \\log\_2 N$ para imágenes de $N \\times N$).



**Frecuencia de Transformada (Transform Frequency):**

En cualquier transformada unitaria $\\mathbf{y} = \\mathbf{A}\\mathbf{x}$, el concepto de frecuencia se generaliza contando el número de cruces por cero de las filas (funciones base) de la matriz $\\mathbf{A}$.

Si las filas de $\\mathbf{A}$ se ordenan de forma que los cruces por cero aumenten con el número de fila, los coeficientes resultantes $y(k)$ quedan ordenados según su número de onda o frecuencia de transformada creciente (por ejemplo, frecuencia del coseno).

Es un concepto relativo a cada transformada: un término de baja frecuencia en una transformada puede contener armónicos de alta frecuencia de otra.


**La Transformada Óptima (The Optimum Transform):**

Para evaluar el rendimiento en filtrado y compresión de datos según el criterio del error cuadrático medio (MSE).



#### 6.3 Las propiedades de las transformaciones unitarias - conservación de la energía y rotación

**Conservación de la Energía (1D):**

Una transformación unitaria $\\mathbf{v} = \\mathbf{A}\\mathbf{u}$ conserva intacta la energía de la señal (la norma al cuadrado o longitud del vector $\\mathbf{u}$ en el espacio $N$-dimensional):



&#x09;	$$\\Vert{}\\mathbf{v}\\Vert{}^2 = \\Vert{}\\mathbf{u}\\Vert{}^2$$



dado que $\\mathbf{A}^{\*T}\\mathbf{A} = \\mathbf{I}$)



**Interpretación Geométrica (Rotación):**

Al no alterar la longitud del vector, toda transformación unitaria equivale simplemente a una rotación del vector $\\mathbf{u}$ en el espacio $N$-dimensional (o de forma equivalente, una rotación de los ejes coordenados base, donde las componentes de $\\mathbf{v}$ son las proyecciones de $\\mathbf{u}$ sobre la nueva base).



**Conservación de la Energía en 2D:**



La misma propiedad se cumple para imágenes bidimensionales: la suma de las energías de los píxeles en el dominio espacial es exactamente igual a la suma de las energías de los coeficientes en el dominio de la transformada:



&#x09;	$$\\sum\_{m,n=0}^{N-1} \\vert{}u(m,n)\\vert{}^2 = \\sum\_{k,l=0}^{N-1} \\vert{}v(k,l)\\vert{}^2$$



#### 6.4. Transformación del Coseno $\\mathbf{C}$

Los elementos $c(k,n)$ de la matriz de la Transformada Discreta del Coseno (DCT) de dimensión $N \\times N$ se definen como:



$$c(k,n) = \\begin{cases} \\dfrac{1}{\\sqrt{N}}, \& k = 0, \\quad 0 \\le n \\le N-1 \\\\ \\sqrt{\\dfrac{2}{N}} \\cos \\dfrac{\\pi(2n+1)k}{2N}, \& 1 \\le k \\le N-1, \\quad 0 \\le n \\le N-1 \\end{cases}$$



**Definición de la DCT Unidimensional:** Para una secuencia $\\{u(n), \\, 0 \\le n \\le N-1\\}$, su DCT directa en 1D se calcula mediante:



&#x09;	$$v(k) = \\alpha(k) \\sum\_{n=0}^{N-1} u(n) \\cos\\left\[\\frac{\\pi(2n+1)k}{2N}\\right], \\quad 0 \\le k \\le N-1$$



Donde el factor de normalización $\\alpha(k)$ garantiza la ortonormalidad de la transformación:



$$\\alpha(0) \\triangleq \\sqrt{\\frac{1}{N}}, \\qquad \\alpha(k) \\triangleq \\sqrt{\\frac{2}{N}} \\quad \\text{para } 1 \\le k \\le N-1$$



**Transformada Inversa del Coseno en 1D**

La secuencia original $u(n)$ se reconstruye a partir de los coeficientes $v(k)$ mediante:



&#x09;	$$u(n) = \\sum\_{k=0}^{N-1} \\alpha(k) v(k) \\cos\\left\[\\frac{\\pi(2n+1)k}{2N}\\right], \\quad 0 \\le n \\le N-1$$



**Par de Transformadas del Coseno en 2D:**

Dado que la matriz $\\mathbf{C}$ es real y ortogonal ($\\mathbf{A} = \\mathbf{A}^\* = \\mathbf{C}$), el par de transformadas 2D se obtiene sustituyendo directamente $\\mathbf{C}$ en las ecuaciones matriciales separables.



**Directa: $\\mathbf{V} = \\mathbf{C}\\mathbf{U}\\mathbf{C}^T$**

**Inversa: $\\mathbf{U} = \\mathbf{C}^T\\mathbf{V}\\mathbf{C}$**



**Propiedades:**

* La matriz de la transformada del coseno $\\mathbf{C}$ está compuesta únicamente por números reales ($\\mathbf{C} = \\mathbf{C}^\*$) y es ortogonal, por lo que su inversa es simplemente su matriz transpuesta (sin necesidad de conjugación compleja):



&#x09;	$$\\mathbf{C} = \\mathbf{C}^\* \\implies \\mathbf{C}^{-1} = \\mathbf{C}^T$$



* La DCT no es la parte real de la Transformada Discreta de Fourier (DFT) unitaria $\\mathbf{F}$. Sin embargo, la DCT de una secuencia sí está directamente relacionada con la DFT de su extensión simétrica.



**Reconstrucción de la Imagen (DCT Inversa):**

Para recuperar la imagen aproximada en el dominio espacial, se aplica la transformada inversa aprovechando que la matriz es ortogonal ($\\mathbf{C}^{-1} = \\mathbf{C}^T$):



&#x09;	$$\\mathbf{U} = \\mathbf{C}^T\\mathbf{V}\\mathbf{C}$$



Esto equivale a reconstruir la imagen como una suma ponderada únicamente de las imágenes base más representativas ($\\mathbf{U} = \\sum v(k,l)\\mathbf{A}\_{k,l}^\*$).



**Criterio del Error Cuadrático Medio (MSE) y Relación Señal a Ruido (SNR)**

Para evaluar cuantitativamente la fidelidad de una imagen reconstruida $u'(m,n)$ respecto a la imagen original $u(m,n)$ de dimensiones $M \\times N$ tras un proceso de compresión, se emplea el criterio del error cuadrático medio.



Desde un enfoque estadístico, el error cuadrático medio (mean square error, $\\sigma\_{ms}^2$) se define mediante la esperanza matemática de la diferencia cuadrática entre ambas imágenes:



&#x09;	$$\\sigma\_{ms}^2 \\triangleq E\\left\[\\vert{}u(m,n) - u'(m,n)\\vert{}^2\\right]$$



Asimismo, el error cuadrático medio promedio ($\\sigma\_a^2$) sobre toda la extensión espacial de la imagen está dado por:



&#x09;	$$\\sigma\_a^2 = \\frac{1}{MN} \\sum\_{m=1}^{M} \\sum\_{n=1}^{N} E\\left\[\\vert{}u(m,n) - u'(m,n)\\vert{}^2\\right]$$



En aplicaciones prácticas, cuando no se dispone del conjunto estadístico completo, se utiliza el error cuadrático promedio (average least squares error, $\\sigma\_{ls}^2$) como estimador directo:



&#x09;	$$\\sigma\_{ls}^2 \\triangleq \\frac{1}{MN} \\sum\_{m=1}^{M} \\sum\_{n=1}^{N} \\vert{}u(m,n) - u'(m,n)\\vert{}^2$$



En virtud de la propiedad de conservación de energía de las transformadas unitarias (teorema de Parseval), el error cuadrático evaluado en el dominio espacial es exactamente equivalente al error cuadrático calculado entre los coeficientes en el dominio de la transformada.



En el ámbito del procesamiento y codificación de imágenes, el error cuadrático medio ($\\sigma\_e^2$, donde $\\sigma\_e$ representa a $\\sigma\_a$, $\\sigma\_{ms}$ o $\\sigma\_{ls}$) suele expresarse en decibelios ($\\text{dB}$) mediante la Relación Señal a Ruido ($\\text{SNR}$).



&#x09;	$$\\text{SNR} = 10 \\log\_{10} \\frac{\\sigma^2}{\\sigma\_e^2}$$



donde $\\sigma^2$ corresponde a la varianza de la imagen original de referencia.



De manera alternativa, en aplicaciones específicas de compresión y codificación de imágenes se emplea la Relación Señal a Ruido Pico ($\\text{SNR}'$ o $\\text{PSNR}$), la cual normaliza el error respecto al cuadrado del valor pico a pico máximo de la imagen de referencia:



&#x09;	$$\\text{SNR}' = 10 \\log\_{10} \\frac{(\\text{valor pico a pico de la imagen de referencia})^2}{\\sigma\_e^2}$$





## **7. Referencias y Bibliografía**

**\[1]** V. Britanak, P. C. Yip, and K. R. Rao, Discrete Cosine and Sine Transforms: General Properties, Fast Algorithms and Integer Approximations. Oxford, U.K.: Academic Press, 2007

**\[2]** A. K. Jain, Fundamentals of Digital Image Processing. Englewood Cliffs, NJ, USA: Prentice-Hall, 1989.









