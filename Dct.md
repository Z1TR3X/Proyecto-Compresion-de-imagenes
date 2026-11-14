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



## **7. Referencias y Bibliografía**

**\[1]** V. Britanak, P. C. Yip, and K. R. Rao, Discrete Cosine and Sine Transforms: General Properties, Fast Algorithms and Integer Approximations. Oxford, U.K.: Academic Press, 2007









