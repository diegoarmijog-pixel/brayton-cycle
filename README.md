# Simulador de Ciclo Brayton con Oxicombustión y Captura de CO₂

Simulador termodinámico avanzado para análisis de ciclos Brayton con oxicombustión, CO₂ supercrítico y captura de carbono.

---

## Características Principales

### 🔥 Oxicombustión Pura
- **Sin emisiones de N₂**: Oxígeno puro de ASU (Air Separation Unit)
- **Productos 100% CO₂ + H₂O**: Facilita captura de carbono
- **Recirculación de CO₂**: Control térmico y dilución de llama

### ⚙️ Múltiples Combustibles
- **Gas Natural**: 95% CH₄, 3% C₂H₆, 1% C₃H₈, 1% N₂
- **Metano Puro**: 100% CH₄
- **Gas de Síntesis**: Relaciones H₂/CO variables (1:1, 2:1, etc.)
- **Hidrógeno**: 100% H₂
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
- **Exportación CSV**: Resultados descargables para análisis

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
     ├→ C7: H₂O líquida (captura)
     └→ C8: CO₂ seco
         ↓
        [División]
         ├→ C9: CO₂ a captura (1-f_recir)
         │   ↓
         │  [Compresor CO₂ Captura]  ← W_comp_captura
         │   ↓
         │  C11: CO₂ almacenamiento (P_storage)
         │
         └→ C12: CO₂ a recirculación (f_recir)
             ↓
            [Compresor CO₂ Recirc]  ← W_comp_recir
             ↓
            C8: CO₂ comprimido
             ↓
            [Recuperador - Lado Frío]
             ↓
            C10: CO₂ precalentado
             └──→ [vuelve a combustión]
```

### Balance Energético

```
W_neto = W_turbina - W_comp1 - W_comp_recir - W_comp_captura - W_ASU

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
| `matplotlib` | latest | Gráficos base |
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
- **Presión Almacenamiento CO₂**: 50-300 bar (default: 150 bar)

### 3. Ejecutar Simulación

Presionar botón **"🚀 SIMULAR"**

El simulador ejecuta:
1. Cálculo estequiométrico de combustión
2. Balance energético (determina T_combustión)
3. Expansión en turbina (trabajo generado)
4. Proceso de separación de agua
5. Sistema de recirculación y captura
6. Balances energéticos totales

### 4. Resultados

#### Pestaña 1: Diagrama del Proceso
- Diagrama visual del ciclo
- **Balance Energético**:
  - Trabajo turbina (MW)
  - Trabajos de compresión (MW)
  - Trabajo ASU (MW)
  - Trabajo neto (MW)
  - Calor de combustión (MW)
  - **Eficiencia térmica (%)**
- **Flujos de CO₂**:
  - CO₂ total generado
  - CO₂ recirculado
  - CO₂ capturado
  - Fracción de recirculación efectiva

#### Pestaña 2: Propiedades Termodinámicas
Tabla completa con 13 corrientes:
- Nombre de corriente
- Temperatura (°C, K)
- Presión (bar, MPa)
- Entalpía específica (kJ/mol)
- Entropía específica (J/(mol·K))
- Flujo molar (mol/s)
- Flujo másico (kg/s)
- Composición molar (%)

Botón de **descarga CSV** para análisis externo.

#### Pestaña 3: Diagrama P-H
Gráfico interactivo Presión vs. Entalpía con trayectoria del ciclo.

#### Pestaña 4: Diagrama T-S
Gráfico interactivo Temperatura vs. Entropía con trayectoria del ciclo.

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

| Combustible | LHV (J/mol) | Fuente |
|-------------|-------------|--------|
| Gas Natural | 802000 | NIST |
| Metano | 802000 | NIST |
| Syngas (1:1) | 250000 | Calculado conservador |
| Hidrógeno | 241800 | NIST |
| Propano | 2043000 | NIST |
| Etanol | 1277000 | NIST |

---

## Documentación Técnica

### Archivos Disponibles

```
simulador definitivo/
│
├── simulador_brayton.py                      # Código principal (3508 líneas)
├── README.md                                 # Este archivo
├── requirements.txt                          # Dependencias Python
├── diagrama_brayton_corregido.png           # Diagrama del proceso
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
        T, P: Temperatura y presión
        composicion: Dict[str, float] (fracciones molares)
        flujo_molar: float (mol/s)
        h, s, rho: Propiedades calculadas

    Métodos:
        calcular_propiedades(metodo_mezcla=None)
            - None → PR con kij (default)
            - "fugacidad" → HEOS + fugacidad (combustión)
    """

class Combustible:
    """
    Define características de combustibles.

    Atributos:
        tipo: str
        composicion: Dict[str, float]
        LHV: float (J/mol)

    Métodos:
        calcular_productos_combustion(n_O2) → Dict[str, float]
    """

class SimuladorBrayton:
    """
    Motor principal de simulación del ciclo completo.

    Atributos:
        corrientes: List[Corriente] (13 corrientes C0-C12)
        params: Dict (parámetros del ciclo)
        resultados: Dict (W, Q, η)

    Métodos:
        ejecutar_simulacion() → bool
        calcular_temperatura_combustion() → float
        calcular_trabajo_turbina() → float
        calcular_compresores() → Dict[str, float]
        optimizar_parametros() → Dict
    """
```

### Flujo de Ejecución

```
1. Inicialización Streamlit
   ↓
2. Captura de parámetros de usuario
   ↓
3. Crear instancia SimuladorBrayton
   ↓
4. ejecutar_simulacion()
   ├─ Corriente C0: Combustible ambiente
   ├─ Corriente C1: Compresión combustible
   ├─ Corriente C2: O₂ de ASU
   ├─ Balance combustión → T_combustión (iterativo)
   ├─ Corriente C3: Productos combustión
   ├─ Corriente C4: Expansión turbina
   ├─ Corriente C5: Después de recuperador (caliente)
   ├─ Corriente C6: Antes de separador
   ├─ Corriente C7: H₂O líquida
   ├─ Corriente C8: CO₂ seco
   ├─ División flujo
   ├─ Corriente C9: CO₂ a captura
   ├─ Corriente C11: CO₂ almacenamiento
   ├─ Corriente C12: CO₂ a recirculación
   ├─ Corriente C8: CO₂ comprimido (recirc)
   └─ Corriente C10: CO₂ precalentado
   ↓
5. Calcular trabajos y calores
   ↓
6. Calcular eficiencia térmica
   ↓
7. Generar resultados y visualizaciones
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

## Optimización y Mejora Continua

### Funcionalidad de Optimización

El simulador incluye un **optimizador automático** que busca maximizar eficiencia térmica ajustando:

- Presión de combustión (50-250 bar)
- Fracción de recirculación CO₂ (60-95%)
- Efectividad del recuperador (50-95%)

**Algoritmo**: Nelder-Mead simplex (scipy.optimize)

**Restricciones**:
- T_combustión < 2000°C
- Trabajo neto > 0
- Balance energético cerrado

### Casos de Validación

El código incluye validación contra:

1. **Caso CO₂ puro supercrítico** (500 K, 100 bar)
   - Comparación HEOS vs. datos NIST
   - Error h < 0.05%, s < 0.1%

2. **Caso mezcla CO₂-H₂O** (600 K, 150 bar, x_CO₂=0.7)
   - Comparación PR+kij vs. datos experimentales
   - Error ρ < 3%, h < 2%

3. **Ciclo Brayton simple** (solo CH₄, sin recirculación)
   - Comparación con ciclo Brayton clásico
   - Eficiencia coherente con literatura

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

### Versión 3.0 (Noviembre 2024)
- ✅ Implementación de kij experimental para CO₂-H₂O (0.1896)
- ✅ Eliminación de ecuaciones SRK e ideal (código limpiado)
- ✅ Documentación consolidada actualizada
- ✅ Método dual: HEOS para puros, PR para mezclas
- ✅ Correcciones de exceso complementarias a kij
- ✅ Optimizador de parámetros del ciclo
- ✅ Validación con datos experimentales

### Versión 2.0 (Noviembre 2024)
- Método HEOS + fugacidad para combustión
- Correcciones CO₂-H₂O de exceso
- Sistema de 13 corrientes completo
- Diagramas P-H y T-S interactivos

### Versión 1.0 (Octubre 2024)
- Versión inicial del simulador
- Interfaz Streamlit básica
- Cálculos termodinámicos fundamentales

---

**🚀 Desarrollado para el avance de tecnologías de captura de carbono y ciclos de potencia supercríticos**
