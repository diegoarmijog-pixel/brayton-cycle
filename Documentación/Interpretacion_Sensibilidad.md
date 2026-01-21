# Interpretación Termodinámica del Análisis de Sensibilidad

## Relación x_H₂O(C4) vs f_recirculación

- **↑ Recirculación** → Más CO₂ recirculado en C3 → **Diluye productos de combustión** → **↓ x_H₂O** en C4
- **↓ Recirculación** → Menos dilución por CO₂ → Mayor fracción de productos (H₂O + CO₂ fresco) → **↑ x_H₂O** en C4

## ¿Por qué η_Global aumenta hasta un punto y luego baja?

### Efecto 1 - Trabajo del Compresor CO₂ (negativo)
- ↓ x_H₂O significa ↑ recirculación → ↑ W_comp_CO₂ (más CO₂ a comprimir desde baja P hasta alta P)
- Esto **REDUCE** η_Global (más trabajo parásito consumido)

### Efecto 2 - Potencia de la Turbina (positivo)
- ↑ Recirculación → ↑ flujo másico total en turbina (más moles de gas) → ↑ W_turbina
- Esto **AUMENTA** η_Global (más potencia generada)

### Efecto 3 - Recuperación de Calor (positivo)
- ↑ Recirculación → Mayor flujo de CO₂ precalentado en recuperador → Reduce necesidad de combustible
- Esto **AUMENTA** η_Global (mejor aprovechamiento térmico de gases de escape)

### Efecto 4 - Temperatura de Combustión (restricción operacional)
- **↓ x_H₂O** (↑ recirculación) → ↓ T_combustión (dilución térmica por CO₂)
- **↑ x_H₂O** (↓ recirculación) → ↑ T_combustión (menos dilución, más concentración de calor)
- **RESTRICCIÓN FÍSICA:** T_combustión debe ser < 1800-2000°C para evitar daño a materiales (turbina, cámara)
- Por lo tanto, **aunque η_Global sea mayor a la derecha del gráfico (alta x_H₂O), esas condiciones son INOPERABLES** debido a T_combustión excesiva

## Conclusión
El **óptimo real** no es el máximo de η_Global, sino el **máximo dentro de la zona operable** (donde T_combustión < 1800°C). Revisar la columna 'T_combustión (°C)' en la tabla de datos numéricos para identificar el rango operacional. Típicamente, esto ocurre en la zona de **baja x_H₂O** (< 30-40%), correspondiente a **alta recirculación** (> 85-90%).