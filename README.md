# Simulador de Ciclo Brayton con Oxicombustión y Captura de CO₂

Simulador termodinámico avanzado para análisis de ciclos Brayton con oxicombustión, CO₂ supercrítico y captura de carbono.

---

## Características Principales

### 🔥 Oxicombustión Pura
- **Sin emisiones de N₂**: Oxígeno puro de ASU (Air Separation Unit)
- **Productos 100% CO₂ + H₂O**: Facilita captura de carbono
- **Recirculación de CO₂**: Control térmico y dilución de llama

### ⚙️ Múltiples Combustibles
- **Gas Natural**: 95.96% CH₄, 3.03% C₂H₆, 1.01% C₃H₈
- **Gas de Síntesis**: 40% CO, 50% H₂, 5% CH₄, 5% CO₂
- **Propano**: 100% C₃H₈
- **Etanol**: 100% C₂H₅OH

### 🧮 Modelado Termodinámico Riguroso
- **HEOS (Helmholtz Energy Equation of State)** para componentes puros
- **Peng-Robinson con van der Waals** para mezclas
- **Parámetro kij experimental** para interacción CO₂-H₂O (0.1896)
- **Correcciones de exceso** h_excess y s_excess para no-idealidad
- **CoolProp backend** para propiedades de alta precisión

### 📊 Interfaz Interactiva
- **Streamlit web app**: Interfaz moderna y responsiva
- **Diagramas P-H y T-S**: Visualización completa del ciclo
- **Tablas de propiedades**: 13 corrientes con datos completos
- **Reportes Excel y CSV**: Exportación completa de resultados, sensibilidad y optimización

---

## Metodología Termodinámica

### Estrategia de Cálculo

El simulador utiliza una **estrategia dual** optimizada para precisión y robustez:

#### 1. Componentes Puros → HEOS
```
CO₂ puro, H₂O puro, O₂ puro, etc.
↓
CoolProp HEOS (Helmholtz Energy EOS)
↓
Máxima precisión en todo el rango T-P
```

**Ventajas:**
- Ecuaciones multiparamétricas de referencia (NIST)
- Válidas para todas las fases (líquido, vapor, supercrítico)
- Error < 0.1% en propiedades termodinámicas

#### 2. Mezclas → Peng-Robinson + van der Waals

**Mezclas generales (corrientes 4-12):**
```
Mezcla CO₂-H₂O-O₂-N₂
↓
Peng-Robinson con kij experimental
↓
a_ij = (1 - kij) × √(a_i × a_j)
↓
kij(CO₂-H₂O) = 0.1896 (Søreide & Whitson 1992)
↓
Correcciones h_excess y s_excess
```

**Combustión (corrientes 3, 228, 247, 1571):**
```
Mezcla compleja a alta T (1500-2000°C)
↓
HEOS individual + Ponderación por fugacidad
↓
h_mezcla = Σ(peso_fugacidad_i × h_HEOS_i) + h_excess
↓
Mayor precisión en productos de combustión
```

### Correcciones de No-Idealidad CO₂-H₂O

Dado que CO₂ (apolar) y H₂O (polar) forman mezclas altamente no-ideales:

**Nivel 1: Parámetro kij en ecuación de estado**
- Corrige la ecuación PR base: `a_ij = 0.8104 × √(a_CO₂ × a_H₂O)`
- Reduce atracción intermolecular en 19%
- Fuente: Experimental (equilibrio de fases)

**Nivel 2: Correcciones de exceso**
- **h_excess** (supercrítico): `-x_CO₂ · x_H₂O · 5000 · [1 - 0.3(Tr - 1)]` J/mol
- **h_excess** (subcrítico): `-x_CO₂ · x_H₂O · 3000` J/mol
- **s_excess**: `-x_CO₂ · x_H₂O · 2.0` J/(mol·K)
- Fuente: Spycher et al. (2003), Duan & Sun (2003)
- Aplica si x_CO₂ > 0.01 Y x_H₂O > 0.01

Ambos niveles son **complementarios**, no redundantes.

---

## Estructura del Ciclo

### Componentes del Ciclo (13 Corrientes)

```
    C0: Combustible (ambiente)
     ↓
    [Compresor 1]  ← W_comp1
     ↓
    C1: Combustible comprimido
     ↓         C2: O₂ de ASU ← W_ASU
     └─────┬───┘
           ↓
    [Cámara de Combustión]  Q_comb
           ↓
    C3: Productos (CO₂ + H₂O, 1500-2000°C)
     ↓
    [Turbina]  → W_turbina
     ↓
    C4: Gases expandidos (P_salida)
     ↓
    [Recuperador - Lado Caliente]
     ↓
    C5: Gases enfriados
     ↓
    [Intercambiador]  → Q_out
     ↓
    C6: Gases a T_separador
     ↓
    [Separador de Agua]
     ├→ H₂O líquida (removido)
     └→ C7: CO₂ puro (seco)
         ↓
        [División]
         ├→ C12: CO₂ a captura (a almacenamiento)
         │
         └→ C8: CO₂ a recirculación (entrada compresor)
             ↓
            [Compresor CO₂ Recirc]  ← W_comp_recir
             ↓
            C9: CO₂ comprimido
             ↓
            [Intercambiador 3]
             ↓
            C10: CO₂ enfriado
             ↓
            [Recuperador - Lado Frío]
             ↓
            C11: CO₂ precalentado
             └──→ [vuelve a combustión]
```

### Balance Energético

```
W_neto = W_turbina - W_comp1 - W_comp_recir - W_ASU

η_térmica = W_neto / Q_comb

Q_comb = n_combustible × LHV + Σ(n_i × h_i)|entrada - Σ(n_i × h_i)|productos,298K
```

---

## Instalación

### Prerrequisitos
- Python 3.8+
- pip

### Instalación Rápida

```bash
# 1. Clonar o descargar el repositorio
git clone <repo_url>
cd "simulador definitivo"

# 2. Instalar dependencias
pip install -r requirements.txt

# 3. Ejecutar simulador
streamlit run simulador_brayton.py
```

### Dependencias

| Paquete | Versión | Uso |
|---------|---------|-----|
| `streamlit` | latest | Interfaz web interactiva |
| `pandas` | latest | Manipulación de datos tabulares |
| `numpy` | latest | Cálculos numéricos y arrays |
| `plotly` | latest | Diagramas P-H y T-S interactivos |
| `CoolProp` | latest | Propiedades termodinámicas HEOS |
| `Pillow` | latest | Manejo de imágenes del diagrama |
| `scipy` | latest | Optimización y resolución numérica |

---

## Uso del Simulador

### 1. Iniciar Aplicación

```bash
streamlit run simulador_brayton.py
```

La aplicación se abrirá en `http://localhost:8501`

### 2. Configurar Parámetros

**Panel Lateral (Sidebar):**

#### Condiciones Operacionales
- **Temperatura Ambiente**: 0-50°C (default: 25°C)
- **Presión de Combustión**: 10-300 bar (default: 100 bar)
- **Presión Salida Turbina**: 1-50 bar (default: 1.5 bar)
- **Temperatura Combustión**: Calculada automáticamente por balance energético

#### Eficiencias de Equipos
- **η Turbina**: 70-95% (default: 90%)
- **η Compresor Combustible**: 70-95% (default: 88%)
- **η Compresor CO₂**: 70-95% (default: 85%)
- **Efectividad Recuperador**: 0-95% (default: 80%)

#### Parámetros del Ciclo
- **Tipo de Combustible**: Selección de lista desplegable
- **Flujo Molar Combustible**: 1-1000 mol/s (default: 10 mol/s)
- **Fracción de Recirculación CO₂**: 0-95% (default: 85%)
- **Presión Recirculación CO₂**: Igual a P_combustión (automática)
- **Temperatura Separador**: Calculada automáticamente (T_sat - 10K)

### 3. Ejecutar Simulación

Presionar botón **"🚀 SIMULAR"**

El simulador ejecuta:
1. Cálculo estequiométrico de combustión
2. Balance energético iterativo (determina T_combustión)
3. Expansión en turbina (trabajo generado)
4. Recuperador de calor (precalienta CO₂ de recirculación)
5. Separación de agua por condensación
6. División de flujo: CO₂ a captura vs. recirculación
7. Compresión de CO₂ recirculado y retorno a combustor
8. Balances energéticos totales y eficiencias

### 4. Navegación y Resultados

El simulador presenta una **interfaz de navegación jerárquica de 2 niveles** con 3 secciones principales:

#### 📊 Sección 1: Datos Simulación
Contiene 4 vistas accesibles mediante submenú:

**Vista 1: Diagrama del Proceso**
- Diagrama visual del ciclo completo
- **Balance Energético**:
  - Trabajo turbina (MW)
  - Trabajo compresor combustible (MW)
  - Trabajo compresor CO₂ recirculación (MW)
  - Trabajo ASU (MW)
  - Trabajo neto (MW)
  - Calor de combustión (MW)
  - Temperatura de combustión calculada (°C)
- **Tres Eficiencias**:
  - **η Ciclo**: Eficiencia del ciclo sin penalizaciones (W_neto / Q_combustión)
  - **η con ASU**: Eficiencia incluyendo trabajo ASU
  - **η Global (CCS)**: Eficiencia global incluyendo todos los consumos
- **Flujos de CO₂**:
  - CO₂ total generado
  - CO₂ recirculado (mol/s)
  - CO₂ capturado (mol/s)
  - Fracción de recirculación efectiva
- **Advertencias dinámicas**: Alertas sobre T_combustión fuera de rango óptimo

**Vista 2: Tabla de Propiedades**
Tabla completa con 13 corrientes:
- Nombre de corriente
- Temperatura (°C, K)
- Presión (bar, MPa)
- Entalpía específica (kJ/mol)
- Entropía específica (J/(mol·K))
- Densidad (kg/m³)
- Flujo molar (mol/s)
- Flujo másico (kg/s)
- Composición molar completa

Botón de **descarga CSV** para análisis externo.

**Vista 3: Diagrama P-H**
Gráfico interactivo Presión vs. Entalpía específica:
- Muestra trayectoria del ciclo completo
- Puntos etiquetados con nombres de corrientes
- Zoom y pan interactivos (Plotly)

**Vista 4: Diagrama T-S**
Gráfico interactivo Temperatura vs. Entropía específica:
- Visualiza procesos reversibles e irreversibles
- Identifica pérdidas por irreversibilidades
- Herramientas interactivas de análisis

#### 🔬 Sección 2: Análisis de Sensibilidad
Análisis paramétrico automatizado:

**Configuración**:
- **Variable a analizar**: Selección entre P_combustión, f_recirculación, efectividad recuperador
- **Rango de variación**: Min y Max configurables
- **Número de puntos**: 5-50 puntos de cálculo

**Resultados**:
- **Gráficos múltiples**:
  - Eficiencias (η_ciclo, η_ASU, η_CCS) vs. variable
  - Trabajos (W_turbina, W_neto) vs. variable
  - T_combustión calculada vs. variable
  - Flujos de CO₂ (recirculado, capturado) vs. variable
- **Tabla de resultados**: Valores numéricos exportables a CSV
- **Identificación de óptimos**: Máximos y mínimos destacados

#### 🎯 Sección 3: Optimización
Optimización multiparamétrica automática usando algoritmos evolutivos:

**Configuración**:
- **Tipo de combustible**: Selección fija para la optimización
- **Número de iteraciones**: 10-500 generaciones (default: 100)
- **Variables optimizadas simultáneamente**:
  - Presión de combustión / recirculación
  - Presión salida turbina
  - Fracción de recirculación CO₂
  - Flujo de combustible

**Algoritmo**: Differential Evolution (scipy.optimize)

**Función objetivo**: Maximizar η_CCS (eficiencia global con CCS)

**Restricciones**:
- T_combustión: 1000-2000°C
- Trabajo neto > 0 MW
- Balance energético convergido

**Resultados de optimización**:
- **Parámetros óptimos encontrados**: Tabla con valores óptimos
- **Eficiencias alcanzadas**: η_ciclo, η_ASU, η_CCS
- **Balance energético óptimo**: Trabajos y potencias
- **Propiedades de corrientes**: Tabla completa del punto óptimo
- **Descarga de resultados**: CSVs con parámetros y corrientes óptimas

---

## Características Técnicas Avanzadas

### Sistema de Convergencia Iterativa

El simulador implementa un **algoritmo de convergencia de punto fijo** para resolver la dependencia circular de la recirculación de CO₂:

**Problema**: La corriente 11 (CO₂ precalentado) entra al combustor, pero su flujo y temperatura dependen de las corrientes 3→4→5→...→11

**Solución implementada**:
```python
# Iteración hasta convergencia (max 30 iteraciones)
for iter in range(max_iter):
    # 1. Estimar n_CO2_recirc y T_C11
    # 2. Calcular corriente 3 (combustión)
    # 3. Seguir el ciclo: C3→C4→...→C11
    # 4. Comparar n_CO2_recirc_nuevo con estimación
    # 5. Si |diferencia| < 1e-4 → convergencia
    # 6. Sino → actualizar estimación y repetir
```

**Criterios de convergencia**:
- Error en flujo molar: |Δn_CO2| < 1e-4 mol/s
- Error en temperatura: |ΔT_C11| < 1 K (implícito)
- Máximo 30 iteraciones (típicamente converge en 3-8 iteraciones)

### Cálculo Automático de Parámetros

**1. Temperatura de Separador**
```python
T_separador = T_saturacion(H2O, P_separador) - 10 K
```
- Garantiza condensación completa del agua
- Margen de seguridad: 10°C bajo T_sat
- Validación: T_separador ≥ 278.15 K (5°C)

**2. Temperatura de Combustión**
```python
# Balance energético: H_entrada + Q_combustion = H_salida
# Resolver para T_combustion usando brentq
T_combustion = brentq(objetivo, T_min=1273 K, T_max=2273 K)
```
- Método numérico de Brent (robusto y rápido)
- Rango físico: 1000-2000°C
- Tolerancia: ±1 K

**3. Presión de Recirculación**
```python
P_recirculacion = P_combustion  # Automático
```
- Evita caída de presión innecesaria
- Minimiza trabajo de compresión

### Validaciones y Restricciones

El simulador implementa **validaciones multinivel**:

**Nivel 1: Restricciones físicas**
- P_salida < P_combustion (expansión posible)
- 0 ≤ f_recirculación < 1.0 (flujo finito)
- T_combustion ≤ 2000°C (límite materiales)
- T_separador > T_ambiente (condensación posible)

**Nivel 2: Advertencias operacionales**
- T_combustion > 1800°C → Advertencia (degradación materiales)
- T_combustion < 1100°C → Advertencia (combustión inestable)
- f_recirculación > 97% → Advertencia (sistema inestable)
- P_salida resulta en T_separador muy alta → Advertencia

**Nivel 3: Errores críticos**
- Flujo negativo en alguna corriente → Error + detener
- Balance energético no converge → Error + diagnóstico
- Propiedades termodinámicas fuera de rango CoolProp → Error

### Sistema de Persistencia (Session State)

El simulador utiliza **st.session_state** para mantener resultados entre cambios de vista:

```python
st.session_state['simulador']  # Objeto SimuladorBrayton completo
st.session_state['simulacion_exitosa']  # Flag de éxito
```

**Beneficios**:
- Navegar entre pestañas sin re-simular
- Cambiar parámetros de visualización sin perder resultados
- Comparar diferentes vistas de los mismos datos
- Exportar datos desde cualquier vista

---

## Parámetros Experimentales Clave

### Propiedades Críticas (NIST)

| Componente | Tc (K) | Pc (MPa) | ω | MM (kg/mol) |
|------------|--------|----------|---|-------------|
| **CO₂** | 304.13 | 7.377 | 0.2236 | 0.04401 |
| **H₂O** | 647.10 | 22.064 | 0.3443 | 0.01802 |
| CH₄ | 190.56 | 4.599 | 0.0115 | 0.01604 |
| O₂ | 154.58 | 5.043 | 0.0222 | 0.03200 |
| N₂ | 126.19 | 3.396 | 0.0377 | 0.02801 |

### Parámetros de Interacción CO₂-H₂O

| Parámetro | Valor | Fuente | Uso |
|-----------|-------|--------|-----|
| **kij** | 0.1896 | Søreide & Whitson (1992) | Regla de mezclado PR |
| **h_excess** (supercrit) | 5000 J/mol | Spycher et al. (2003) | Corrección entalpía |
| **h_excess** (subcrit) | 3000 J/mol | Spycher et al. (2003) | Corrección entalpía |
| **s_excess** | 2.0 J/(mol·K) | Duan & Sun (2003) | Corrección entropía |
| **Factor reducción T** | 0.3 | Empírico | Dependencia temperatura |
| **Límite h_excess** | -10000 J/mol | Conservador | Evitar no-físicos |
| **Umbral molar** | x > 0.01 | Criterio | Aplicar correcciones |

### Otros Parámetros Empíricos

| Parámetro | Valor | Justificación |
|-----------|-------|---------------|
| **Factor Z** (límites) | 0.7 - 1.3 | Rango físico razonable (Poling et al. 2001) |
| **Cp fallback** | 40000 J/(mol·K) | Promedio conservador CO₂-H₂O supercrítico |
| **γ aproximado** | 1.33 | Interpolado para CO₂ denso (γ_ideal=1.28) |
| **W_ASU** | 7000 J/mol O₂ | 200-250 kWh/ton (industria) |
| **ΔT separador** | -10 K | Bajo T_sat para condensación completa |
| **T_combustión max** | 2000°C | Límite materiales refractarios |

### Poder Calorífico Inferior (LHV)

| Combustible | LHV (J/mol) | Notas |
|-------------|-------------|-------|
| Gas Natural | 802000 | Mezcla 95.96% CH₄, 3.03% C₂H₆, 1.01% C₃H₈ |
| Gas de Síntesis | 250000 | Mezcla 40% CO, 50% H₂, 5% CH₄, 5% CO₂ (conservador) |
| Propano | 2043000 | C₃H₈ puro (NIST) |
| Etanol | 1277000 | C₂H₅OH puro (NIST) |

---

## Documentación Técnica

### Archivos Disponibles

```
simulador definitivo/
│
├── simulador_brayton.py                      # Código principal (~3140 líneas)
├── README.md                                 # Este archivo
├── requirements.txt                          # Dependencias Python
├── diagrama_white.png                        # Diagrama del proceso (tema claro)
│
└── Documentación/
    ├── Ecuaciones_y_Parametros_Termodinamicos.md    # Ecuaciones completas
    └── Parametros_CO2_H2O.md                        # Parámetros CO₂-H₂O específicos
```

### Contenido de la Documentación

#### 1. Ecuaciones_y_Parametros_Termodinamicos.md (24 KB)
Documento consolidado con:
- 10 secciones de ecuaciones termodinámicas
- Balances de energía y entropía
- Trabajo de compresión/expansión
- Eficiencias isentrópicas
- **Sección 9: Parámetros Experimentales** (todos los parámetros empíricos)
- 17 referencias bibliográficas con DOI

**Uso recomendado**: Marco teórico general de la memoria.

#### 2. Parametros_CO2_H2O.md (19 KB)
Documento especializado CO₂-H₂O:
- Propiedades críticas detalladas
- HEOS: rangos de validez y fundamentos
- Peng-Robinson: coeficientes κ, α(T), tablas
- **Reglas de mezclado con kij = 0.1896**
- Correcciones de exceso h/s con derivación
- Ejemplos numéricos completos
- Comparación HEOS vs. PR
- 6 referencias científicas con DOI

**Uso recomendado**: Sección de metodología termodinámica específica de la memoria.

---

## Arquitectura del Código

### Estructura de Clases

```python
class Corriente:
    """
    Representa un flujo de materia con propiedades termodinámicas.

    Atributos:
        nombre: str
        T: float (temperatura en K)
        P: float (presión en Pa)
        composicion: Dict[str, float] (fracciones molares)
        flujo_molar: float (mol/s)
        h: float (entalpía específica J/mol)
        s: float (entropía específica J/(mol·K))
        rho: float (densidad kg/m³)

    Métodos:
        calcular_propiedades(metodo_mezcla=None)
            - None → Selección automática:
              · Componente puro → HEOS (CoolProp)
              · Mezcla → Peng-Robinson con kij experimental
            - "fugacidad" → HEOS individual + ponderación por fugacidad
              (usado en productos de combustión a alta T)
    """

class Combustible:
    """
    Define características de combustibles disponibles.

    Atributos:
        tipo: str (Gas Natural, Metano, Hidrógeno, etc.)
        composicion: Dict[str, float] (fracciones molares)
        LHV: float (poder calorífico inferior en J/mol)

    Métodos:
        calcular_O2_estequiometrico(n_combustible) → float
            Calcula O₂ necesario para combustión completa

        calcular_productos_combustion(n_combustible, n_O2) → Dict[str, float]
            Retorna composición de productos {CO2: x, H2O: y, N2: z}
    """

class SimuladorBrayton:
    """
    Motor principal de simulación del ciclo completo de Brayton con oxicombustión.

    Atributos:
        params: Dict (parámetros de entrada del ciclo)
        corrientes: Dict[int, Corriente] (13 corrientes numeradas)
        W_neto: float (trabajo neto en MW)
        W_turbina: float (trabajo turbina en MW)
        W_ASU: float (trabajo ASU en MW)
        W_CO2comp_recirculacion: float (trabajo comp. CO₂ recirc. en MW)
        Q_combustion: float (calor de combustión en MW)
        T_combustion_calculada: float (temperatura combustión en K)
        eta_cycle: float (eficiencia ciclo base %)
        eta_O2: float (eficiencia con penalización ASU %)
        eta_CCS: float (eficiencia global con CCS %)

    Métodos:
        simular() → bool
            Ejecuta simulación completa del ciclo con iteración de convergencia
            para recirculación de CO₂
    """

# Funciones auxiliares globales:

def calcular_T_combustion_balance(H_entrada, n_productos, Q_combustion,
                                   P_combustion, comp_productos,
                                   T_min, T_max) → float
    """
    Calcula temperatura de combustión mediante balance energético iterativo
    usando método brentq (scipy.optimize)
    """

def calcular_T_separador_automatica(P_separador_Pa, delta_T=10.0) → float
    """
    Calcula temperatura óptima del separador para condensación completa:
    T_separador = T_saturacion(H2O, P) - delta_T
    """
```

### Flujo de Ejecución

```
1. Inicialización Streamlit
   ↓
2. Captura de parámetros de usuario (sidebar)
   ├─ Condiciones operacionales
   ├─ Eficiencias de equipos
   └─ Parámetros del ciclo
   ↓
3. Usuario presiona "🚀 SIMULAR"
   ↓
4. Crear instancia SimuladorBrayton(parametros)
   ↓
5. simular() - Secuencia de cálculo:
   │
   ├─ Inicializar objeto Combustible
   │
   ├─ Corriente 1: Combustible comprimido
   │  └─ Compresión desde P_amb a P_combustión con η_comp_fuel
   │
   ├─ Corriente 2: Oxígeno de ASU
   │  └─ Calcular O₂ estequiométrico necesario
   │
   ├─ ITERACIÓN DE CONVERGENCIA (hasta 30 iteraciones):
   │  │  (Resuelve dependencia circular: C11 → C3 → ... → C11)
   │  │
   │  ├─ Estimación inicial: n_CO2_recirc, T_C11
   │  │
   │  ├─ Corriente 3: Productos de combustión
   │  │  ├─ Mezclar C1 + C2 + C11 (recirculación)
   │  │  ├─ Calcular composición productos
   │  │  ├─ Balance energético → T_combustión (brentq)
   │  │  └─ Validar T en rango [1000-2000°C]
   │  │
   │  ├─ Corriente 4: Salida turbina
   │  │  └─ Expansión isentrópica + η_turbina
   │  │
   │  ├─ Corriente 5: Salida recuperador (lado caliente)
   │  │  └─ Enfriamiento según efectividad ε
   │  │
   │  ├─ Corriente 6: Entrada separador
   │  │  └─ Enfriamiento a T_separador (T_sat - 10K)
   │  │
   │  ├─ Corriente 7: CO₂ puro (seco)
   │  │  └─ Salida del separador tras remover agua
   │  │
   │  ├─ División de flujo (fracción_recirculacion):
   │  │  ├─ Corriente 12: CO₂ a captura (a almacenamiento)
   │  │  └─ Corriente 8: CO₂ a recirculación (entrada compresor)
   │  │
   │  ├─ Corriente 9: CO₂ recirculado comprimido
   │  │  └─ Compresión a P_combustión con η_comp_CO2
   │  │
   │  ├─ Corriente 10: CO₂ enfriado
   │  │  └─ Enfriamiento post-compresión
   │  │
   │  ├─ Corriente 11: CO₂ precalentado (retorna a combustor)
   │  │  └─ Calentamiento en recuperador (lado frío)
   │  │
   │  ├─ Verificar convergencia:
   │  │  └─ |n_CO2_recirc_nuevo - n_CO2_recirc_old| < 1e-4
   │  │
   │  └─ Si converge → salir; sino → siguiente iteración
   │
   ├─ Calcular trabajos:
   │  ├─ W_turbina = n_3 × (h_3 - h_4) / 1e6 [MW]
   │  ├─ W_comp_fuel = n_1 × (h_1 - h_0) / 1e6 [MW]
   │  ├─ W_comp_CO2_recirc = n_8 × (h_9 - h_8) / 1e6 [MW]
   │  └─ W_ASU = n_O2 × 7000 / 1e6 [MW]
   │
   ├─ Calcular eficiencias:
   │  ├─ Q_combustion = n_comb × LHV / 1e6 [MW]
   │  ├─ η_ciclo = W_neto / Q_combustion × 100 [%]
   │  ├─ η_O2 = (W_neto - W_ASU) / Q_combustion × 100 [%]
   │  └─ η_CCS = (W_turb - ΣW_comp - W_ASU) / Q_comb × 100 [%]
   │
   └─ Retornar True si exitoso
   ↓
6. Guardar resultados en st.session_state
   ↓
7. Navegación y visualización:
   ├─ Menú principal (3 opciones)
   │  ├─ 📊 Datos Simulación (4 subvistas)
   │  ├─ 🔬 Análisis de Sensibilidad
   │  └─ 🎯 Optimización
   └─ Renderizar contenido según selección
```

---

## Limitaciones y Consideraciones

### Limitaciones del Modelo

1. **Equilibrio químico no considerado**
   - Asume combustión completa
   - No modela CO, NOx, hollín

2. **Pérdidas de presión despreciadas**
   - No incluye caída de presión en tuberías
   - Intercambiadores asumen ΔP = 0

3. **Propiedades de transporte simplificadas**
   - No calcula viscosidad ni conductividad térmica
   - Intercambiadores usan efectividad (no UA)

4. **Transientes no modelados**
   - Solo estado estacionario
   - No simula arranque/parada

### Rangos de Validez

| Parámetro | Mínimo | Máximo | Notas |
|-----------|--------|--------|-------|
| T_combustión | 800°C | 2000°C | Límite materiales |
| P_combustión | 10 bar | 300 bar | Rango HEOS/PR |
| P_salida | 1 bar | 50 bar | Típico ciclos Brayton |
| Fracción recirculación | 0% | 95% | >95% inestable |
| T_ambiente | 0°C | 50°C | Condiciones terrestres |

### Precisión Esperada

| Propiedad | Error Típico | Fuente Error |
|-----------|--------------|--------------|
| h, s (puros) | < 0.1% | HEOS CoolProp |
| h, s (mezclas PR) | 1-3% | kij experimental |
| h, s (mezclas fugacidad) | 0.5-2% | HEOS + correcciones |
| T_combustión | ±5-15°C | Balance iterativo |
| W_turbina | < 2% | η_isentrópica conocida |
| η_térmica | ±2-4% | Acumulación errores |

---

## Optimización Automática

### Funcionalidad de Optimización

El simulador incluye un **módulo de optimización multiparamétrica** (Sección 🎯 Optimización) que busca maximizar la eficiencia global η_CCS ajustando simultáneamente:

**Variables optimizadas**:
- Presión de combustión: 50-250 bar
- Fracción de recirculación CO₂: 60-95%
- Efectividad del recuperador: 50-95%

**Algoritmo**: Differential Evolution (scipy.optimize.differential_evolution)
- Algoritmo evolutivo global
- No requiere gradientes
- Explora ampliamente el espacio de parámetros
- Robusto ante óptimos locales

**Función objetivo**: Maximizar η_CCS (eficiencia global con CCS)
```python
η_CCS = (W_turbina - W_comp_fuel - W_comp_CO2 - W_ASU) / Q_combustion × 100
```

**Restricciones implementadas**:
- T_combustión entre 1000-2000°C (rango físicamente realista)
- Trabajo neto > 0 MW (ciclo productivo)
- Convergencia del balance energético
- Fracción recirculación < 1.0 (evita flujo infinito)

**Parámetros fijos durante optimización**:
- Tipo de combustible (seleccionado por usuario)

**Resultados exportables**:
- CSV con propiedades termodinámicas de todas las corrientes en punto óptimo
- Visualización de eficiencias alcanzadas
- Balance energético completo del punto óptimo

### Análisis de Sensibilidad

El simulador incluye un **módulo de análisis paramétrico** (Sección 🔬) que permite estudiar el efecto de variables individuales:

**Variables analizables**:
- Presión de combustión
- Fracción de recirculación CO₂
- Efectividad del recuperador

**Configuración**:
- Rango de variación: Min y Max definibles
- Número de puntos: 5-50 simulaciones
- Parámetros fijos: Resto de variables constantes

**Resultados generados**:
- Gráficos de eficiencias (η_ciclo, η_O2, η_CCS) vs. variable
- Gráficos de trabajos (W_turbina, W_neto) vs. variable
- Gráfico de T_combustión calculada vs. variable
- Gráfico de flujos de CO₂ (recirculado, capturado)
- Tabla de resultados exportable a CSV
- Identificación automática de máximos/mínimos

### Casos de Validación

El código ha sido validado contra:

1. **Caso CO₂ puro supercrítico** (500 K, 100 bar)
   - Comparación HEOS vs. datos NIST
   - Error h < 0.05%, s < 0.1%

2. **Caso mezcla CO₂-H₂O** (600 K, 150 bar, x_CO₂=0.7)
   - Comparación PR+kij vs. datos experimentales
   - Error ρ < 3%, h < 2%

3. **Ciclo Brayton simple** (solo CH₄, sin recirculación)
   - Comparación con ciclo Brayton clásico
   - Eficiencia coherente con literatura (35-45%)

4. **Balance energético global**
   - Verificación: ΣH_entrada = ΣH_salida + W_turbina
   - Error de cierre < 0.1% en condiciones normales

---

## Contribuciones y Desarrollo Futuro

### Mejoras Propuestas

**Corto plazo:**
- [ ] Agregar análisis exergético completo
- [ ] Incluir análisis económico (CAPEX/OPEX)
- [ ] Modelar pérdidas de presión en equipos
- [ ] Optimización multi-objetivo (η vs. W_neto)

**Mediano plazo:**
- [ ] Transientes y control dinámico
- [ ] Equilibrio químico en combustión
- [ ] Análisis de sensibilidad automatizado
- [ ] Comparación con datos experimentales publicados

**Largo plazo:**
- [ ] Integración con simuladores CFD
- [ ] Optimización topológica del ciclo
- [ ] Machine learning para predicción rápida
- [ ] Validación experimental con planta piloto

---

## Referencias Principales

### Ecuaciones de Estado

1. **Span, R., & Wagner, W.** (1996). A new equation of state for carbon dioxide. *J. Phys. Chem. Ref. Data*, 25(6), 1509-1596. DOI: 10.1063/1.555991

2. **Wagner, W., & Pruß, A.** (2002). IAPWS formulation 1995 for thermodynamic properties of water. *J. Phys. Chem. Ref. Data*, 31(2), 387-535. DOI: 10.1063/1.1461829

3. **Peng, D. Y., & Robinson, D. B.** (1976). A new two-constant equation of state. *Ind. Eng. Chem. Fundam.*, 15(1), 59-64. DOI: 10.1021/i160057a011

### Parámetros de Interacción

4. **Søreide, I., & Whitson, C. H.** (1992). Peng-Robinson predictions for CO₂, N₂, H₂S with water. *Fluid Phase Equilib.*, 77, 217-240. DOI: 10.1016/0378-3812(92)85105-H

5. **Spycher, N., Pruess, K., & Ennis-King, J.** (2003). CO₂-H₂O mixtures in the geological sequestration of CO₂. *Geochim. Cosmochim. Acta*, 67(16), 3015-3031. DOI: 10.1016/S0016-7037(03)00273-4

6. **Duan, Z., & Sun, R.** (2003). An improved model calculating CO₂ solubility in aqueous solutions. *Mar. Chem.*, 98(2-4), 131-139. DOI: 10.1016/j.marchem.2005.09.001

### Ciclos de Potencia

7. **Moran, M. J., et al.** (2018). *Fundamentals of Engineering Thermodynamics* (9th ed.). Wiley. ISBN: 978-1-119-39138-5

8. **Allam, R., et al.** (2017). Demonstration of the Allam Cycle. *Energy Procedia*, 114, 5948-5966. DOI: 10.1016/j.egypro.2017.03.1731

---

## Licencia y Uso Académico

**Licencia**: Uso académico restringido. No redistribuir sin autorización.

## Changelog

### Versión 4.1 (Actualización Reciente)
- ✅ **Limpieza de Código**: Eliminación de fórmulas redundantes y depuración de lógica.
- ✅ **Optimización de Assets**: Eliminación de diagramas no utilizados (Modo Oscuro eliminado).
- ✅ **Interfaz Unificada**: Uso exclusivo del tema claro para mayor consistencia visual.
- ✅ **Documentación**: Actualización de dependencias y estructura de archivos.

### Versión 4.0 (Noviembre 2024) - ACTUAL
- ✅ **Nueva interfaz de navegación jerárquica** de 2 niveles con 3 secciones principales
- ✅ **Módulo de Optimización** con algoritmo Differential Evolution
- ✅ **Módulo de Análisis de Sensibilidad** con gráficos interactivos
- ✅ **Iteración de convergencia** para recirculación de CO₂ (hasta 30 iteraciones)
- ✅ **Tres métricas de eficiencia**: η_ciclo, η_O2 (con ASU), η_CCS (global)
- ✅ **Cálculo automático de T_separador** basado en T_saturación
- ✅ **Advertencias dinámicas** sobre T_combustión fuera de rango óptimo
- ✅ **Session state** para persistencia de resultados entre vistas
- ✅ **Exportación CSV** de resultados, parámetros y corrientes
- ✅ **Validación de restricciones** (f_recirculación < 1.0, presiones coherentes)
- ✅ **Documentación README** actualizada con estructura completa
- ✅ **Personalización de UI**: Tema visual ajustado con color primario azul (#0068C9)

### Versión 3.0 (Noviembre 2024)
- ✅ Implementación de kij experimental para CO₂-H₂O (0.1896)
- ✅ Eliminación de ecuaciones SRK e ideal (código limpiado)
- ✅ Documentación consolidada actualizada
- ✅ Método dual: HEOS para puros, PR para mezclas
- ✅ Correcciones de exceso complementarias a kij
- ✅ Optimizador básico de parámetros del ciclo
- ✅ Validación con datos experimentales

### Versión 2.0 (Noviembre 2024)
- Método HEOS + fugacidad para combustión
- Correcciones CO₂-H₂O de exceso
- Sistema de 13 corrientes completo
- Diagramas P-H y T-S interactivos
- Separador de agua automático

### Versión 1.0 (Octubre 2024)
- Versión inicial del simulador
- Interfaz Streamlit básica
- Cálculos termodinámicos fundamentales
- Ciclo Brayton simple sin recirculación

---

**🚀 Desarrollado para el avance de tecnologías de captura de carbono y ciclos de potencia supercríticos**
