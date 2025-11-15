# Ecuaciones Termodinámicas y Parámetros Fisicoquímicos
## Simulador de Ciclo Brayton con Oxicombustión y Captura de CO₂

Este documento describe todas las ecuaciones termodinámicas, parámetros fisicoquímicos y constantes utilizadas en el simulador de ciclo Brayton abierto con oxicombustión.

---

## Tabla de Contenidos

1. [Constantes Universales](#1-constantes-universales)
2. [Propiedades Críticas de Componentes](#2-propiedades-críticas-de-componentes)
3. [Sistema de Ecuaciones de Estado](#3-sistema-de-ecuaciones-de-estado)
4. [Ecuaciones de Estado para Mezclas](#4-ecuaciones-de-estado-para-mezclas)
5. [Correcciones de No-Idealidad](#5-correcciones-de-no-idealidad)
6. [Balance de Masa y Energía](#6-balance-de-masa-y-energía)
7. [Trabajo de Compresores y Turbina](#7-trabajo-de-compresores-y-turbina)
8. [Eficiencias del Ciclo](#8-eficiencias-del-ciclo)
9. [Parámetros Experimentales y Ajustes](#9-parámetros-experimentales-y-ajustes)
10. [Referencias Bibliográficas](#10-referencias-bibliográficas)

---

## 1. Constantes Universales

### 1.1 Constantes Fundamentales

| Parámetro | Símbolo | Valor | Unidad | Referencia |
|-----------|---------|-------|--------|------------|
| Constante universal de los gases | R | 8.314 | J/(mol·K) | CODATA 2018 |
| Temperatura de referencia estándar | T₀ | 298.15 | K | IUPAC |
| Presión de referencia estándar | P₀ | 101325 | Pa | IUPAC |

**Referencias:**
- CODATA 2018: DOI: 10.1103/RevModPhys.93.025010
- IUPAC Green Book (2007)

---

## 2. Propiedades Críticas de Componentes

### 2.1 Componentes Principales del Ciclo

| Componente | Fórmula | Tc (K) | Pc (MPa) | ω | MM (kg/mol) | Referencia |
|------------|---------|--------|----------|---|-------------|------------|
| **CO₂** | CO₂ | 304.13 | 7.377 | 0.2236 | 0.04401 | Span & Wagner (1996) |
| **H₂O** | H₂O | 647.10 | 22.064 | 0.3443 | 0.01802 | Wagner & Pruß (2002) |
| **O₂** | O₂ | 154.58 | 5.043 | 0.0222 | 0.03200 | Schmidt & Wagner (1985) |

### 2.2 Componentes de Combustible

| Componente | Fórmula | Tc (K) | Pc (MPa) | ω | MM (kg/mol) | Referencia |
|------------|---------|--------|----------|---|-------------|------------|
| **CH₄** | CH₄ | 190.56 | 4.599 | 0.0115 | 0.01604 | Setzmann & Wagner (1991) |
| **C₂H₆** | C₂H₆ | 305.32 | 4.872 | 0.0995 | 0.03007 | Bücker & Wagner (2006) |
| **C₃H₈** | C₃H₈ | 369.83 | 4.248 | 0.1523 | 0.04410 | Lemmon et al. (2009) |
| **C₂H₅OH** | C₂H₅OH | 513.92 | 6.148 | 0.6450 | 0.04607 | Dillon & Penoncello (2004) |
| **CO** | CO | 132.86 | 3.494 | 0.0497 | 0.02801 | Lemmon & Span (2006) |
| **H₂** | H₂ | 33.19 | 1.313 | -0.2160 | 0.00202 | Leachman et al. (2009) |

**Nota:** Factor acéntrico negativo del H₂ se debe a efectos cuánticos.

### 2.3 Poder Calorífico Inferior (LHV)

| Combustible | LHV (kJ/mol) | LHV (MJ/kg) | Referencia |
|-------------|--------------|-------------|------------|
| Gas Natural (CH₄ equiv.) | 802 | ~50.0 | NIST WebBook |
| Gas de Síntesis* | 250 | ~10-12 | Calculado |
| Propano (C₃H₈) | 2,043 | ~46.4 | NIST WebBook |
| Etanol (C₂H₅OH) | 1,277 | ~27.7 | NIST WebBook |

*Gas de síntesis: 40% CO + 50% H₂ + 5% CH₄ + 5% CO₂ (base molar)

---

## 3. Sistema de Ecuaciones de Estado

### 3.1 Jerarquía de Selección

El simulador selecciona automáticamente el método termodinámico según:

```
┌─────────────────────────────────────┐
│  ¿Es componente puro o mezcla?      │
└─────────────┬───────────────────────┘
              │
        ┌─────┴─────┐
        │           │
      PURO        MEZCLA
        │           │
        ▼           ▼
     HEOS      Peng-Robinson
  (CoolProp)   + van der Waals
```

### 3.2 HEOS para Componentes Puros

**Ecuación:** Helmholtz Energy Equation of State (multiparamétrica)

**Forma general:**
$$a(T, \rho) = a^{ideal}(T, \rho) + a^{res}(T, \rho)$$

Donde:
- $a = A/(RT)$: Energía libre de Helmholtz reducida
- $a^{ideal}$: Contribución del gas ideal
- $a^{res}$: Contribución residual (interacciones moleculares)

**Propiedades derivadas:**
$$P = \rho RT\left(1 + \delta \frac{\partial a^{res}}{\partial \delta}\right)$$
$$h = RT\left[1 + \tau\left(\frac{\partial a^{ideal}}{\partial \tau} + \frac{\partial a^{res}}{\partial \tau}\right) + \delta\frac{\partial a^{res}}{\partial \delta}\right]$$
$$s = R\left[\tau\left(\frac{\partial a^{ideal}}{\partial \tau} + \frac{\partial a^{res}}{\partial \tau}\right) - a^{ideal} - a^{res}\right]$$

Donde:
- $\delta = \rho/\rho_c$: Densidad reducida
- $\tau = T_c/T$: Temperatura reducida inversa

**Rango de validez:**
- CO₂: 216.6 K - 1100 K, hasta 800 MPa
- H₂O: 273.16 K - 2000 K, hasta 1000 MPa
- Otros: Según base de datos CoolProp

**Implementación:** `CoolProp.AbstractState("HEOS", "CO2")`

---

## 4. Ecuaciones de Estado para Mezclas

### 4.1 Método Peng-Robinson (PR)

**Ecuación cúbica:**
$$P = \frac{RT}{V_m - b} - \frac{a(T)}{V_m^2 + 2bV_m - b^2}$$

**Parámetros para componente puro i:**

$$a_i = 0.45724 \frac{R^2 T_{c,i}^2}{P_{c,i}} \alpha_i(T)$$

$$b_i = 0.07780 \frac{RT_{c,i}}{P_{c,i}}$$

$$\alpha_i = [1 + \kappa_i(1-\sqrt{T_r,i})]^2$$

$$\kappa_i = 0.37464 + 1.54226\omega_i - 0.26992\omega_i^2$$

**Regla de mezclado de van der Waals:**

$$a_{mix} = \sum_i \sum_j x_i x_j \sqrt{a_i a_j}(1 - k_{ij})$$

$$b_{mix} = \sum_i x_i b_i$$

**Parámetro de interacción binaria:** $k_{ij} = 0$ (mezcla ideal, excepto CO₂-H₂O)

**Factor de compresibilidad:**

Resolver ecuación cúbica:
$$Z^3 - (1-B)Z^2 + (A-3B^2-2B)Z - (AB-B^2-B^3) = 0$$

Donde:
$$A = \frac{a_{mix}P}{(RT)^2}, \quad B = \frac{b_{mix}P}{RT}$$

**Seleccionar raíz máxima real** (fase vapor)

### 4.2 Funciones de Desviación (Departure Functions)

**Entalpía residual:**
$$h^{dep} = RT\left[Z - 1 - \frac{A}{2\sqrt{2}B}\ln\left(\frac{Z+(1+\sqrt{2})B}{Z+(1-\sqrt{2})B}\right)\right]$$

**Entropía residual:**
$$s^{dep} = R\left[\ln(Z-B) - \frac{A}{2\sqrt{2}B}\ln\left(\frac{Z+(1+\sqrt{2})B}{Z+(1-\sqrt{2})B}\right)\right]$$

**Propiedades finales:**
$$h_{mix} = h^{ideal} + h^{dep}$$
$$s_{mix} = s^{ideal} + s^{dep} - R\sum_i x_i \ln(x_i)$$

Donde:
- $h^{ideal} = \sum_i x_i h_i(T, P)$ calculado con HEOS por componente
- $s^{ideal} = \sum_i x_i s_i(T, P)$ calculado con HEOS por componente
- $-R\sum_i x_i \ln(x_i)$: Corrección de entropía de mezclado ideal

### 4.3 Densidad de la Mezcla

$$\rho_{molar} = \frac{P}{ZRT}$$

$$\rho_{másica} = \rho_{molar} \cdot MM_{mix}$$

$$MM_{mix} = \sum_i x_i MM_i$$

---

## 5. Correcciones de No-Idealidad

### 5.1 Factor de Compresibilidad Acotado

Para evitar valores no físicos en mezclas con HEOS:

$$Z_{acotado} = \max(0.7, \min(1.3, Z_{calculado}))$$

**Justificación física:**
- $Z < 0.7$: Líquidos altamente comprimidos (no relevante en gases de combustión)
- $Z \approx 1.0$: Gas ideal
- $Z > 1.3$: Poco común para mezclas CO₂-H₂O hasta 300 bar

**Referencia:** Poling et al. (2001) - *The Properties of Gases and Liquids*

### 5.2 Correcciones de Interacción Binaria CO₂-H₂O

Aplicables cuando $x_{CO_2} > 0.01$ y $x_{H_2O} > 0.01$

#### 5.2.1 Entalpía de Exceso

**CO₂ supercrítico** ($T_r > 1.0$, donde $T_r = T/T_{c,CO_2}$):
$$h^{excess} = -x_{CO_2} \cdot x_{H_2O} \cdot 5000 \cdot [1 - 0.3(T_r - 1)] \text{ J/mol}$$

Con límite: $h^{excess} \geq -10000$ J/mol

**CO₂ subcrítico** ($T_r \leq 1.0$):
$$h^{excess} = -x_{CO_2} \cdot x_{H_2O} \cdot 3000 \text{ J/mol}$$

#### 5.2.2 Entropía de Exceso

$$s^{excess} = -x_{CO_2} \cdot x_{H_2O} \cdot 2.0 \text{ J/(mol·K)}$$

**Aplicación:**
$$h_{final} = h_{mix} + h^{excess}$$
$$s_{final} = s_{mix} + s^{excess}$$

**Referencia:** Spycher et al. (2003) - DOI: 10.1016/S0016-7037(03)00273-4

### 5.3 Regla de Mezclado Basada en Fugacidad

Para HEOS con mezclas complejas, se usa peso de fugacidad:

$$w_i = \frac{x_i \cdot f_i}{\sum_j x_j \cdot f_j}$$

Donde $f_i$ es la fugacidad del componente puro i a T, P.

**Entalpía de mezcla con peso de fugacidad:**
$$h_{mix} = \sum_i w_i h_i$$

**Referencia:** Prausnitz et al. (1999) - *Molecular Thermodynamics*

---

## 6. Balance de Masa y Energía

### 6.1 Estequiometría de Combustión

#### Metano (Gas Natural):
$$\text{CH}_4 + 2\text{O}_2 \rightarrow \text{CO}_2 + 2\text{H}_2\text{O}$$

#### Propano:
$$\text{C}_3\text{H}_8 + 5\text{O}_2 \rightarrow 3\text{CO}_2 + 4\text{H}_2\text{O}$$

#### Etanol:
$$\text{C}_2\text{H}_5\text{OH} + 3\text{O}_2 \rightarrow 2\text{CO}_2 + 3\text{H}_2\text{O}$$

#### Gas de Síntesis:
$$\text{CO} + 0.5\text{O}_2 \rightarrow \text{CO}_2$$
$$\text{H}_2 + 0.5\text{O}_2 \rightarrow \text{H}_2\text{O}$$

### 6.2 Coeficientes Estequiométricos

| Combustible | Fórmula | O₂ requerido | CO₂ generado | H₂O generado |
|-------------|---------|--------------|--------------|--------------|
| Metano | CH₄ | 2.0 | 1.0 | 2.0 |
| Etano | C₂H₆ | 3.5 | 2.0 | 3.0 |
| Propano | C₃H₈ | 5.0 | 3.0 | 4.0 |
| Etanol | C₂H₅OH | 3.0 | 2.0 | 3.0 |
| Monóxido de carbono | CO | 0.5 | 1.0 | 0.0 |
| Hidrógeno | H₂ | 0.5 | 0.0 | 1.0 |

**Unidades:** mol/mol de combustible

### 6.3 Balance de Energía en la Cámara de Combustión

**Ecuación general:**
$$\dot{Q}_{in} + \sum_{ent} \dot{n}_i h_i = \sum_{sal} \dot{n}_j h_j$$

**Cálculo de temperatura de combustión:**

Dado $\dot{n}_{fuel}$, $\dot{n}_{O_2}$, $\dot{n}_{CO_2,recirc}$, se resuelve iterativamente:

$$\sum_{productos} \dot{n}_j h_j(T_{comb}, P) = \sum_{reactivos} \dot{n}_i h_i(T_{ent}, P) + \dot{n}_{fuel} \cdot LHV$$

**Método numérico:** Bisección sobre T con tolerancia 0.1 K

**Rango de búsqueda:** 800 K - 3500 K

### 6.4 Cálculo de LHV para Gas de Síntesis

Para mezcla de combustibles:
$$LHV_{mix} = \sum_i x_i \cdot LHV_i$$

**Gas de Síntesis estándar** (40% CO, 50% H₂, 5% CH₄, 5% CO₂):
- LHV_CO = 283 kJ/mol
- LHV_H₂ = 242 kJ/mol
- LHV_CH₄ = 802 kJ/mol
- LHV_CO₂ = 0 kJ/mol

$$LHV_{syngas} = 0.40(283) + 0.50(242) + 0.05(802) + 0.05(0) = 274.3 \text{ kJ/mol}$$

**Valor usado (conservador):** 250 kJ/mol

---

## 7. Trabajo de Compresores y Turbina

### 7.1 Turbina (C3 → C4)

**Trabajo ideal (isentrópico):**
$$W_{turb,ideal} = \dot{n}(h_3 - h_{4s})$$

Donde $h_{4s}$ se calcula con $s_{4s} = s_3$ y $P_4 = P_{salida,turbina}$

**Trabajo real:**
$$W_{turb,real} = \eta_{turb} \cdot W_{turb,ideal}$$

**Eficiencia isentrópica típica:** $\eta_{turb} = 0.90$ (referencia: Dixon & Hall, 2014)

### 7.2 Compresor de Combustible (C0 → C1)

**Trabajo ideal:**
$$W_{comp,fuel,ideal} = \dot{n}_{fuel}(h_{1s} - h_0)$$

Donde $h_{1s}$ se calcula con $s_{1s} = s_0$ y $P_1 = P_{combustion}$

**Trabajo real:**
$$W_{comp,fuel,real} = \frac{W_{comp,fuel,ideal}}{\eta_{comp,fuel}}$$

**Eficiencia isentrópica típica:** $\eta_{comp,fuel} = 0.88$ (referencia: Dixon & Hall, 2014)

### 7.3 Compresor de CO₂ Recirculación (C8 → C9)

**Trabajo ideal:**
$$W_{comp,CO_2,ideal} = \dot{n}_{CO_2,recirc}(h_{9s} - h_8)$$

Donde $h_{9s}$ se calcula con $s_{9s} = s_8$ y $P_9 = P_{combustion}$

**Trabajo real:**
$$W_{comp,CO_2,real} = \frac{W_{comp,CO_2,ideal}}{\eta_{comp,CO_2}}$$

**Eficiencia isentrópica típica:** $\eta_{comp,CO_2} = 0.85$ (referencia: Moullec, 2013)

**Nota:** El CO₂ se comprime en condiciones cercanas al punto crítico, lo que reduce la eficiencia respecto a compresores convencionales.

### 7.4 Trabajo del ASU (Air Separation Unit)

**Consumo específico energético:**
$$e_{ASU} = 7000 \text{ J/mol O}_2$$

Equivalente a ~200-250 kWh/ton O₂

**Trabajo total:**
$$W_{ASU} = \dot{n}_{O_2} \cdot e_{ASU}$$

**Referencia:** Valores industriales típicos (Praxair, Air Liquide)

---

## 8. Eficiencias del Ciclo

### 8.1 Eficiencia del Ciclo Básico

Considera solo turbina y compresor de CO₂ recirculación:

$$\eta_{cycle} = \frac{W_{turb} - W_{CO_2,comp,recirc}}{Q_{in}} \times 100\%$$

Donde:
$$Q_{in} = \dot{n}_{fuel} \cdot LHV$$

### 8.2 Eficiencia con ASU

Incluye penalización por producción de oxígeno puro:

$$\eta_{O_2} = \frac{W_{turb} - W_{CO_2,comp,recirc} - W_{ASU}}{Q_{in}} \times 100\%$$

### 8.3 Eficiencia Global (CCS)

Incluye todo el consumo parásito del sistema de captura:

$$\eta_{Global} = \eta_{CCS} = \frac{W_{turb} - W_{CO_2,comp,recirc} - W_{ASU} - W_{comp,fuel}}{Q_{in}} \times 100\%$$

**Trabajo neto del sistema:**
$$W_{neto} = W_{turb} - W_{CO_2,comp,recirc} - W_{ASU} - W_{comp,fuel}$$

**Nota:** Esta es la eficiencia utilizada en la optimización del ciclo.

---

## 9. Parámetros Experimentales y Ajustes

Esta sección contiene parámetros determinados empíricamente o ajustados para el simulador.

### 9.1 Correcciones de Entalpía de Exceso CO₂-H₂O

**Parámetros ajustados:**

| Parámetro | Valor | Unidad | Tipo | Fuente |
|-----------|-------|--------|------|--------|
| Factor base supercrítico | 5000 | J/mol | Experimental | Spycher et al. (2003) |
| Factor temperatura reducida | 0.3 | - | Ajustado | Spycher et al. (2003) |
| Límite inferior | -10000 | J/mol | Conservativo | Spycher et al. (2003) |
| Factor base subcrítico | 3000 | J/mol | Experimental | Spycher et al. (2003) |

**Observación:** Estos valores se determinaron mediante ajuste a datos experimentales de solubilidad y equilibrio de fases CO₂-H₂O en condiciones geológicas (12-100°C, hasta 600 bar).

### 9.2 Parámetro de Interacción Entrópica

| Parámetro | Valor | Unidad | Tipo | Fuente |
|-----------|-------|--------|------|--------|
| s_excess,CO₂-H₂O | 2.0 | J/(mol·K) | Ajustado | Duan & Sun (2003) |

**Observación:** Ajustado a datos de solubilidad de CO₂ en agua salina (273-533 K, 0-2000 bar).

### 9.3 Factor de Compresibilidad Acotado

| Parámetro | Valor | Justificación |
|-----------|-------|---------------|
| Z_min | 0.7 | Evita densidades excesivas en extrapolaciones |
| Z_max | 1.3 | Límite superior físico para gases hasta 300 bar |

**Tipo:** Restricción numérica conservativa

**Fuente:** Poling et al. (2001), observaciones del simulador

### 9.4 Capacidad Calorífica de Respaldo

| Parámetro | Valor | Unidad | Tipo | Justificación |
|-----------|-------|--------|------|---------------|
| Cp_avg | 40000 | J/(mol·K) | Promedio conservador | Valor intermedio para CO₂ y H₂O en rango 300-2000 K |

**Uso:** Fallback cuando CoolProp falla en condiciones extremas.

**Cálculo aproximado:**
$$h \approx 40000(T - 298.15) \text{ J/mol}$$
$$s \approx 40000 \ln(T/298.15) - R\ln(P/P_0) \text{ J/(mol·K)}$$

### 9.5 Consumo Energético del ASU

| Parámetro | Valor | Unidad | Tipo | Fuente |
|-----------|-------|--------|------|--------|
| e_ASU | 7000 | J/mol O₂ | Empírico industrial | Praxair, Air Liquide |
| e_ASU | 200-250 | kWh/ton O₂ | Equivalente | Industria criogénica |

**Observación:** Valor promedio para unidades criogénicas modernas. Puede variar entre 180-300 kWh/ton según pureza y escala.

### 9.6 Eficiencias Isentrópicas de Equipos

| Equipo | Parámetro | Rango típico | Valor por defecto | Tipo | Referencia |
|--------|-----------|--------------|-------------------|------|------------|
| Turbina | η_turb | 0.88-0.92 | 0.90 | Industrial | Dixon & Hall (2014) |
| Comp. combustible | η_comp,fuel | 0.85-0.90 | 0.88 | Industrial | Dixon & Hall (2014) |
| Comp. CO₂ | η_comp,CO₂ | 0.82-0.88 | 0.85 | Experimental | Moullec (2013) |

**Observaciones:**
- Turbina: Valores típicos para turbinas industriales multi-etapa con refrigeración
- Compresor combustible: Compresores centrífugos o axiales convencionales
- Compresor CO₂: Eficiencia reducida por operación cerca del punto crítico

### 9.7 Relación de Calores Específicos

| Sustancia | γ = Cp/Cv | Rango T | Tipo | Fuente |
|-----------|-----------|---------|------|--------|
| CO₂ | 1.33 | 300-500 K | Promedio | CoolProp |

**Uso:** Aproximación rápida en compresión isentrópica:
$$T_{ideal} = T_1 \left(\frac{P_2}{P_1}\right)^{(\gamma-1)/\gamma}$$

**Nota:** En el simulador se usa HEOS completo; este valor es solo para estimaciones iniciales.

### 9.8 Umbral de Aplicación de Correcciones CO₂-H₂O

| Parámetro | Valor | Justificación |
|-----------|-------|---------------|
| x_CO₂,min | 0.01 | Evitar división por cero y efectos numéricos |
| x_H₂O,min | 0.01 | Evitar división por cero y efectos numéricos |

**Tipo:** Restricción numérica

**Aplicación:** Las correcciones de exceso solo se aplican si ambas especies superan 1% molar.

### 9.9 LHV del Gas de Síntesis

| Parámetro | Valor calculado | Valor usado | Tipo | Justificación |
|-----------|-----------------|-------------|------|---------------|
| LHV_syngas | 274.3 kJ/mol | 250 kJ/mol | Conservador | Margen de seguridad para variaciones en composición |

**Cálculo:**
- 40% CO × 283 kJ/mol = 113.2 kJ/mol
- 50% H₂ × 242 kJ/mol = 121.0 kJ/mol
- 5% CH₄ × 802 kJ/mol = 40.1 kJ/mol
- 5% CO₂ × 0 kJ/mol = 0 kJ/mol
- **Total:** 274.3 kJ/mol

**Valor usado:** 250 kJ/mol (8% inferior) para compensar inertes y variabilidad del proceso de gasificación.

### 9.10 Temperatura de Separador

**Cálculo automático:**
$$T_{separador} = T_{sat,H_2O}(P_{salida,turbina}) - 10 \text{ K}$$

**Margen de seguridad:** 10 K por debajo de saturación

**Justificación:** Asegurar condensación completa del H₂O evitando problemas de cavitación.

**Tipo:** Parámetro operacional conservativo

---

## 10. Referencias Bibliográficas

### 10.1 Ecuaciones de Estado

1. **Span, R., & Wagner, W. (1996)**. A new equation of state for carbon dioxide covering the fluid region from the triple‐point temperature to 1100 K at pressures up to 800 MPa. *Journal of Physical and Chemical Reference Data*, 25(6), 1509-1596. DOI: 10.1063/1.555991

2. **Wagner, W., & Pruß, A. (2002)**. The IAPWS formulation 1995 for the thermodynamic properties of ordinary water substance for general and scientific use. *Journal of Physical and Chemical Reference Data*, 31(2), 387-535. DOI: 10.1063/1.1461829

3. **Peng, D. Y., & Robinson, D. B. (1976)**. A new two-constant equation of state. *Industrial & Engineering Chemistry Fundamentals*, 15(1), 59-64. DOI: 10.1021/i160057a011

4. **Bell, I. H., et al. (2014)**. Pure and pseudo-pure fluid thermophysical property evaluation and the open-source thermophysical property library CoolProp. *Industrial & Engineering Chemistry Research*, 53(6), 2498-2508. DOI: 10.1021/ie4033999

### 10.2 Correcciones de No-Idealidad

5. **Spycher, N., Pruess, K., & Ennis-King, J. (2003)**. CO₂-H₂O mixtures in the geological sequestration of CO₂. I. Assessment and calculation of mutual solubilities from 12 to 100°C and up to 600 bar. *Geochimica et Cosmochimica Acta*, 67(16), 3015-3031. DOI: 10.1016/S0016-7037(03)00273-4

6. **Duan, Z., & Sun, R. (2003)**. An improved model calculating CO₂ solubility in pure water and aqueous NaCl solutions from 273 to 533 K and from 0 to 2000 bar. *Chemical Geology*, 193(3-4), 257-271. DOI: 10.1016/S0009-2541(02)00263-2

7. **Prausnitz, J. M., Lichtenthaler, R. N., & de Azevedo, E. G. (1999)**. *Molecular Thermodynamics of Fluid-Phase Equilibria* (3rd ed.). Prentice Hall.

### 10.3 Propiedades de Compuestos Puros

8. **Setzmann, U., & Wagner, W. (1991)**. A new equation of state and tables of thermodynamic properties for methane covering the range from the melting line to 625 K at pressures up to 100 MPa. *Journal of Physical and Chemical Reference Data*, 20(6), 1061-1155. DOI: 10.1063/1.555898

9. **Lemmon, E. W., & Span, R. (2006)**. Short fundamental equations of state for 20 industrial fluids. *Journal of Chemical & Engineering Data*, 51(3), 785-850. DOI: 10.1021/je050186n

10. **Leachman, J. W., et al. (2009)**. Fundamental equations of state for parahydrogen, normal hydrogen, and orthohydrogen. *Journal of Physical and Chemical Reference Data*, 38(3), 721-748. DOI: 10.1063/1.3160306

### 10.4 Bases de Datos y Estándares

11. **NIST Chemistry WebBook**. Linstrom, P. J., & Mallard, W. G. (Eds.). NIST Standard Reference Database Number 69. DOI: 10.18434/T4D303

12. **CODATA 2018**. Tiesinga, E., et al. (2021). CODATA recommended values of the fundamental physical constants: 2018. *Reviews of Modern Physics*, 93(2), 025010. DOI: 10.1103/RevModPhys.93.025010

13. **IUPAC Green Book**. Cohen, E. R., et al. (2007). *Quantities, Units and Symbols in Physical Chemistry* (3rd ed.). RSC Publishing. DOI: 10.1039/9781847557889

### 10.5 Turbomaquinaria y Equipos

14. **Dixon, S. L., & Hall, C. A. (2014)**. *Fluid Mechanics and Thermodynamics of Turbomachinery* (7th ed.). Butterworth-Heinemann. ISBN: 978-0-12-415954-9

15. **Moullec, Y. L. (2013)**. Conceptual study of a high efficiency coal-fired power plant with CO₂ capture using a supercritical CO₂ Brayton cycle. *Energy*, 49, 32-46. DOI: 10.1016/j.energy.2012.10.022

### 10.6 Propiedades de Gases y Mezclas

16. **Poling, B. E., Prausnitz, J. M., & O'Connell, J. P. (2001)**. *The Properties of Gases and Liquids* (5th ed.). McGraw-Hill. ISBN: 0-07-011682-2

17. **Higman, C., & van der Burgt, M. (2008)**. *Gasification* (2nd ed.). Gulf Professional Publishing. ISBN: 978-0-7506-8528-3

---

## Apéndice A: Composiciones de Combustibles

### A.1 Gas Natural (sin N₂)

| Componente | Fracción molar | Base |
|------------|----------------|------|
| CH₄ | 0.9596 | ISO 6976:2016 (normalizado) |
| C₂H₆ | 0.0303 | ISO 6976:2016 (normalizado) |
| C₃H₈ | 0.0101 | ISO 6976:2016 (normalizado) |

**Nota:** Eliminado N₂ (~1%) y redistribuido proporcionalmente para oxicombustión pura.

### A.2 Gas de Síntesis (sin N₂)

| Componente | Fracción molar | Base |
|------------|----------------|------|
| CO | 0.40 | Higman & van der Burgt (2008) |
| H₂ | 0.50 | Higman & van der Burgt (2008) |
| CH₄ | 0.05 | Higman & van der Burgt (2008) |
| CO₂ | 0.05 | Higman & van der Burgt (2008) |

**Nota:** Eliminado N₂ (típico en gasificación con aire) para oxicombustión pura.

---

## Apéndice B: Rangos de Validez

### B.1 Modelos Termodinámicos

| Modelo | T_max (K) | P_max (bar) | Fase | Observaciones |
|--------|-----------|-------------|------|---------------|
| HEOS (CO₂) | 1100 | 8000 | Vapor/Líquido/SC | Máxima precisión |
| HEOS (H₂O) | 2000 | 10000 | Vapor/Líquido/SC | Máxima precisión |
| Peng-Robinson | ~2000 | ~300 | Vapor preferente | Buena precisión en fase gas |
| Fugacidad + HEOS | 2000 | 300 | Mezclas vapor | Robusto para CO₂-H₂O |
| Gas ideal | Sin límite | <10 | Solo vapor | Baja presión únicamente |

### B.2 Condiciones Operacionales Recomendadas

| Parámetro | Rango típico | Restricción física |
|-----------|--------------|-------------------|
| T_combustión | 1200-1800°C | <2000°C (materiales) |
| P_combustión | 10-30 MPa | ≥7.377 MPa (Pc CO₂) |
| P_salida_turbina | 0.1-1.0 MPa | >0.1 MPa (condensación) |
| Fracción recirculación | 85-95% | Control térmico |
| Flujo combustible | 8-25 mol/s | Escala piloto |

---

## Apéndice C: Constantes Numéricas

| Constante | Valor | Uso |
|-----------|-------|-----|
| √2 | 1.41421356 | Departure functions PR |
| ln(10) | 2.302585 | Conversiones logarítmicas |
| Tolerancia bisección | 0.1 K | Cálculo T_combustión |
| Iteraciones máximas | 100 | Convergencia termodinámica |

---

**Documento generado para:** Simulador de Ciclo Brayton con Oxicombustión y Captura de CO₂
**Versión:** 2.0 - Consolidado
**Fecha:** Noviembre 2025
**Autor:** Sistema de Simulación Termodinámica
