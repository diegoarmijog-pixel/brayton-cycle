import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI, AbstractState
from CoolProp import CoolProp as CP
from PIL import Image
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.optimize import differential_evolution


# Configuración de la página
st.set_page_config(page_title="Simulador Ciclo Brayton Oxicombustión", layout="wide")

# Título principal
st.title("Simulador de Ciclo de Brayton Abierto con Oxicombustión")

# ============================================================================
# CLASES Y FUNCIONES AUXILIARES
# ============================================================================

def calcular_T_separador_automatica(P_separador_Pa, delta_T=10.0):
    """
    Calcula la temperatura óptima del separador para garantizar condensación completa del H2O.

    Parámetros:
    -----------
    P_separador_Pa : float
        Presión en el separador en Pa (igual a P_salida_turbina)
    delta_T : float
        Margen de seguridad en °C debajo de la temperatura de saturación (default: 10°C)

    Retorna:
    --------
    T_separador : float
        Temperatura del separador en K
    """
    try:
        # Calcular temperatura de saturación del H2O a la presión del separador
        T_sat_H2O = PropsSI('T', 'P', P_separador_Pa, 'Q', 0, 'Water')  # K

        # T_separador = T_sat - delta_T (para asegurar condensación completa)
        T_separador = T_sat_H2O - delta_T  # K

        # Validación: No puede ser menor que 273.15 K (0°C)
        if T_separador < 273.15:
            T_separador = 273.15 + 5  # Mínimo 5°C

        return T_separador
    except:
        # Si falla (presión muy baja o muy alta), usar valor por defecto
        return 273.15 + 35  # 35°C por defecto


def calcular_T_combustion_balance(H_entrada_combustor, n_total_productos,
                                   Q_combustion, P_combustion,
                                   comp_productos_total, metodo_mezcla=None,
                                   T_min=500+273.15, T_max=3500+273.15):
    """
    Calcula la temperatura de combustión mediante balance energético.

    Balance: H_entrada + Q_combustion = H_salida

    Parámetros:
    -----------
    H_entrada_combustor : float
        Entalpía total de entrada (J/s) = sum(n_i × h_i) para C1, C2, C11
    n_total_productos : float
        Flujo molar total de productos (mol/s)
    Q_combustion : float
        Energía química liberada (J/s) = n_combustible × LHV
    P_combustion : float
        Presión de combustión (Pa)
    comp_productos_total : dict
        Composición molar de productos {"CO2": x, "H2O": y}
    metodo_mezcla : str (opcional, ignorado)
        Mantenido por compatibilidad. El método se selecciona automáticamente:
        - HEOS para componentes puros
        - Peng-Robinson con van der Waals para mezclas
    T_min : float
        Temperatura mínima búsqueda (K) - default 500°C
    T_max : float
        Temperatura máxima búsqueda (K) - default 3500°C (sin límite físico artificial)

    Retorna:
    --------
    T_combustion : float
        Temperatura de combustión calculada (K)
    """
    from scipy.optimize import brentq

    # Entalpía objetivo (J/s)
    H_objetivo = H_entrada_combustor + Q_combustion

    print(f"\n[DEBUG calcular_T_combustion_balance]")
    print(f"  H_entrada_combustor: {H_entrada_combustor/1e6:.4f} MW")
    print(f"  Q_combustion: {Q_combustion/1e6:.4f} MW")
    print(f"  H_objetivo: {H_objetivo/1e6:.4f} MW")
    print(f"  n_total_productos: {n_total_productos:.2f} mol/s")
    print(f"  h_objetivo por mol: {H_objetivo/n_total_productos/1e6:.4f} MJ/mol")
    print(f"  Composición: {comp_productos_total}")
    print(f"  P_combustion: {P_combustion/1e5:.1f} bar")

    # Función objetivo: H_salida(T) - H_objetivo = 0
    def objetivo(T):
        # Crear corriente temporal para calcular h(T, P, comp)
        corriente_temp = Corriente("temp", T=T, P=P_combustion,
                                   composicion=comp_productos_total,
                                   flujo_molar=n_total_productos)
        try:
            # Usar HEOS con mezcla por fugacidad para mayor precisión en combustión
            corriente_temp.calcular_propiedades(metodo_mezcla="fugacidad")
            if corriente_temp.h is None:
                return 1e10  # Error grande
            H_salida = n_total_productos * corriente_temp.h
            return H_salida - H_objetivo
        except:
            return 1e10  # Error grande si falla CoolProp

    # Buscar T que satisfaga el balance
    try:
        # DEBUG: Evaluar en varios puntos
        print(f"\n  Evaluando balance en diferentes temperaturas:")
        T_test = [800, 1000, 1200, 1500, 1800, 2000, 2500, 3000]
        for T_c in T_test:
            T_k = T_c + 273.15
            if T_k >= T_min and T_k <= T_max:
                # Calcular h en este punto
                corr_test = Corriente("test", T=T_k, P=P_combustion,
                                     composicion=comp_productos_total, flujo_molar=n_total_productos)
                corr_test.calcular_propiedades(metodo_mezcla="fugacidad")
                if corr_test.h:
                    H_calc = n_total_productos * corr_test.h
                    error = H_calc - H_objetivo
                    print(f"    T={T_c:4d}°C: h={corr_test.h/1e6:.3f} MJ/mol, H={H_calc/1e6:.2f} MW, error={error/1e6:+.2f} MW")

        T_combustion = brentq(objetivo, T_min, T_max, maxiter=50, xtol=1.0)

        # Validación adicional: verificar que la solución sea físicamente razonable
        error_final = objetivo(T_combustion)
        # print(f"  Solución encontrada: T = {T_combustion-273.15:.1f}°C, error = {error_final/1e6:.4f} MW")

        if abs(error_final) > 1e6:  # Error muy grande
            # El balance no converge bien, usar límites
            if H_objetivo > 0:
                # print(f"  ⚠️ Error grande, retornando T_max")
                return T_max  # Mucha energía disponible
            else:
                # print(f"  ⚠️ Error grande, retornando T_min")
                return T_min  # Poca energía

        return T_combustion
    except ValueError as e:
        # Si no hay solución en el rango, verificar los extremos
        error_min = objetivo(T_min)
        error_max = objetivo(T_max)
        # print(f"  ⚠️ ValueError en brentq: {e}")
        # print(f"  error_min: {error_min/1e6:.2f} MW, error_max: {error_max/1e6:.2f} MW")

        # Si ambos errores son del mismo signo, el rango no contiene la solución
        if (error_min > 0 and error_max > 0):
            # Ambos positivos: H_objetivo es menor que H_salida mínima
            # print(f"  → Combustible insuficiente, retornando T_min")
            return T_min  # Combustible insuficiente
        elif (error_min < 0 and error_max < 0):
            # Ambos negativos: H_objetivo es mayor que H_salida máxima
            # print(f"  → Combustible excesivo, retornando T_max")
            return T_max  # Combustible excesivo
        else:
            # Errores de signos opuestos pero brentq falló
            # print(f"  → Signos opuestos pero brentq falló")
            if abs(error_min) < abs(error_max):
                return T_min
            else:
                return T_max


class Corriente:
    """Clase para representar una corriente en el ciclo"""
    def __init__(self, nombre, T=None, P=None, composicion=None, flujo_molar=None):
        self.nombre = nombre
        self.T = T  # Temperatura en K
        self.P = P  # Presión en Pa
        self.composicion = composicion or {}  # Diccionario {componente: fracción molar}
        self.flujo_molar = flujo_molar  # mol/s
        self.h = None  # Entalpía específica J/mol
        self.s = None  # Entropía específica J/(mol·K)
        self.rho = None  # Densidad kg/m³

    def calcular_propiedades(self, metodo_mezcla=None):
        """
        Calcula las propiedades termodinámicas usando HEOS para puros y PR para mezclas.

        Estrategia de cálculo:
        ----------------------
        - **Componentes puros**: HEOS (Helmholtz Energy Equation of State) vía CoolProp
          - Máxima precisión para todas las fases (líquido, vapor, supercrítico)
          - Cubre todo el rango de T y P relevante para el ciclo

        - **Mezclas**: Peng-Robinson con reglas de van der Waals
          - Parámetro kij experimental para CO₂-H₂O (0.1896)
          - Correcciones de exceso h_excess y s_excess para no-idealidad
          - Robusto para alta temperatura y presión

        Parámetros:
        -----------
        metodo_mezcla : str, optional
            - None (default): Peng-Robinson con van der Waals
            - "fugacidad": HEOS para componentes puros + mezclado por fugacidad
              (más preciso para combustión, usado en cámara de combustión)
        """
        if self.T is None or self.P is None:
            return

        # ====================================================================
        # SELECCIÓN AUTOMÁTICA DE MÉTODO SEGÚN COMPOSICIÓN
        # ====================================================================
        if len(self.composicion) == 1:
            # ================================================================
            # COMPONENTE PURO → Usa HEOS (Helmholtz Equation of State)
            # ================================================================
            # HEOS es el método más preciso para componentes puros
            # Válido para todas las fases (líquido, vapor, supercrítico)
            fluido = list(self.composicion.keys())[0]
            try:
                fluido_coolprop = self._convertir_nombre_fluido(fluido)

                # Usar AbstractState para acceso completo a propiedades
                state = AbstractState("HEOS", fluido_coolprop)
                state.update(CP.PT_INPUTS, self.P, self.T)

                self.h = state.hmass()  # J/kg
                self.s = state.smass()  # J/(kg·K)
                self.rho = state.rhomass()  # kg/m³

                # Convertir a base molar
                MM = state.molar_mass()  # kg/mol
                self.h = self.h * MM  # J/mol
                self.s = self.s * MM  # J/(mol·K)

            except Exception as e:
                # Valores aproximados si falla CoolProp
                self.h = self._h_ideal_gas()
                self.s = self._s_ideal_gas()
                self.rho = self.P / (8.314 * self.T)  # Gas ideal
        else:
            # ================================================================
            # MEZCLA → Método según parámetro
            # ================================================================

            # Si se especifica método de fugacidad, usarlo
            # Este método es más preciso para combustión con HEOS
            if metodo_mezcla == "fugacidad":
                try:
                    self._calcular_mezcla_componentes_puros()
                except Exception as e:
                    # Si falla, usar Peng-Robinson como fallback
                    self._calcular_mezcla_peng_robinson()
            else:
                # Por defecto: Peng-Robinson con van der Waals y kij experimental
                # Robusto para mezclas a alta T y P, especialmente CO2-H2O
                try:
                    self._calcular_mezcla_peng_robinson()
                except Exception as e:
                    # Fallback: método de componentes puros con HEOS
                    self._calcular_mezcla_componentes_puros()

    def _convertir_nombre_fluido(self, nombre):
        """Convierte nombres comunes a nombres de CoolProp"""
        conversion = {
            "CH4": "Methane",
            "C2H6": "Ethane",
            "C3H8": "Propane",
            "H2": "Hydrogen",
            "CO": "CarbonMonoxide",
            "H2O": "Water",
            "CO2": "CarbonDioxide",
            "O2": "Oxygen",
            "N2": "Nitrogen",
            "C2H5OH": "Ethanol"
        }
        return conversion.get(nombre, nombre)

    def _es_mezcla_soportada(self):
        """Verifica si la mezcla es soportada por CoolProp"""
        # CoolProp soporta mezclas predefinidas comunes
        componentes = set(self.composicion.keys())

        # CO2-H2O es una mezcla común en combustión
        if componentes == {"CO2", "H2O"}:
            return True

        # Mezclas con gases comunes
        gases_soportados = {"CO2", "H2O", "N2", "O2", "CH4", "H2", "CO"}
        if componentes.issubset(gases_soportados):
            return True

        return False

    def _calcular_mezcla_componentes_puros(self):
        """
        Calcula propiedades usando HEOS para cada componente puro
        con reglas de mezclado basadas en fugacidad (más robustas)

        ÓPTIMO para CO2 supercrítico + H2O a alta T y P:
        - HEOS individual (muy preciso para cada especie)
        - Reglas de mezclado de fugacidad (no lineales, mejor para sc-CO2)
        - Correcciones de interacción molecular CO2-H2O
        - Válido hasta 2000 K y 300 bar
        """
        R = 8.314  # J/(mol·K)

        # Calcular propiedades individuales y fugacidades
        h_components = {}
        s_components = {}
        Z_components = {}
        MM_components = {}
        fugacity_components = {}

        MM_mezcla = 0

        for componente, fraccion in self.composicion.items():
            if fraccion > 1e-6:
                try:
                    fluido = self._convertir_nombre_fluido(componente)
                    state = AbstractState("HEOS", fluido)
                    state.update(CP.PT_INPUTS, self.P, self.T)

                    h_i = state.hmolar()  # J/mol
                    s_i = state.smolar()  # J/(mol·K)
                    MM_i = state.molar_mass()  # kg/mol

                    # Factor de compresibilidad y fugacidad
                    try:
                        rho_i = state.rhomolar()  # mol/m³
                        Z_i = (self.P / (rho_i * R * self.T)) if rho_i > 0 else 1.0

                        # Fugacidad (indicador de no-idealidad)
                        # f_i = phi_i * P * x_i, donde phi es coef de fugacidad
                        try:
                            fugacity_i = state.fugacity(0) if hasattr(state, 'fugacity') else self.P * Z_i
                        except:
                            fugacity_i = self.P * Z_i
                    except:
                        Z_i = 1.0
                        fugacity_i = self.P

                    h_components[componente] = h_i
                    s_components[componente] = s_i
                    Z_components[componente] = Z_i
                    MM_components[componente] = MM_i
                    fugacity_components[componente] = fugacity_i
                    MM_mezcla += fraccion * MM_i

                except Exception as e:
                    # Fallback conservador
                    try:
                        h_i_mass = PropsSI('H', 'T', self.T, 'P', self.P, fluido)
                        s_i_mass = PropsSI('S', 'T', self.T, 'P', self.P, fluido)
                        MM_i = PropsSI('M', fluido)

                        h_components[componente] = h_i_mass * MM_i
                        s_components[componente] = s_i_mass * MM_i
                        Z_components[componente] = 1.0
                        MM_components[componente] = MM_i
                        fugacity_components[componente] = self.P
                        MM_mezcla += fraccion * MM_i
                    except:
                        pass

        # REGLA DE MEZCLADO BASADA EN FUGACIDAD (no lineal)
        # Para CO2-H2O es superior a reglas lineales en condiciones supercríticas

        self.h = 0
        self.s = 0
        Z_mezcla = 0

        # Ponderación por fugacidad (da más peso a especies con mayor desviación de idealidad)
        total_fugacity_weighted = sum(
            self.composicion[comp] * fugacity_components.get(comp, self.P)
            for comp in self.composicion if self.composicion[comp] > 1e-6 and comp in fugacity_components
        )

        for componente, fraccion in self.composicion.items():
            if fraccion > 1e-6 and componente in h_components:
                # Peso basado en fugacidad (mejor para fases densas)
                f_i = fugacity_components[componente]
                peso_fugacidad = (fraccion * f_i) / total_fugacity_weighted if total_fugacity_weighted > 0 else fraccion

                contribucion_h = peso_fugacidad * h_components[componente]
                self.h += contribucion_h
                self.s += fraccion * s_components[componente]  # Entropía usa fracción molar
                Z_mezcla += fraccion * Z_components[componente]

        # Corrección de entropía de mezclado (termodinámica de mezcla ideal)
        for componente, fraccion in self.composicion.items():
            if fraccion > 1e-6:
                self.s -= R * fraccion * np.log(fraccion)

        # CORRECCIÓN DE INTERACCIÓN MOLECULAR CO2-H2O
        # Parámetros de interacción binaria para mezcla CO2-H2O (literatura)
        if "CO2" in self.composicion and "H2O" in self.composicion:
            x_CO2 = self.composicion.get("CO2", 0)
            x_H2O = self.composicion.get("H2O", 0)

            if x_CO2 > 0.01 and x_H2O > 0.01:
                # Corrección de entalpía de mezcla (no idealidad)
                # Basado en datos experimentales CO2-H2O a alta P y T
                # Interacción débil (sin enlace H) → corrección pequeña pero importante

                Tr_CO2 = self.T / 304.13  # T reducida respecto a Tc del CO2

                # Factor de corrección (empírico, ajustado para sc-CO2 + steam)
                if Tr_CO2 > 1.0:  # CO2 supercrítico
                    # Corrección negativa (mezcla exotérmica débil)
                    h_excess_mix = -x_CO2 * x_H2O * 5000 * (1 - 0.3*(Tr_CO2 - 1))
                    h_excess_mix = max(h_excess_mix, -10000)  # Limitar
                else:
                    h_excess_mix = -x_CO2 * x_H2O * 3000

                self.h += h_excess_mix

                # Corrección de entropía de exceso (pequeña para CO2-H2O)
                s_excess_mix = -x_CO2 * x_H2O * 2.0  # J/(mol·K)
                self.s += s_excess_mix

        # Densidad usando factor Z de la mezcla
        try:
            Z_mezcla = max(0.7, min(1.3, Z_mezcla))  # Rango físico razonable
            self.rho = (self.P * MM_mezcla) / (Z_mezcla * R * self.T)
        except:
            self.rho = 0

    def _calcular_mezcla_peng_robinson(self):
        """
        Calcula propiedades usando ecuación de Peng-Robinson (PR EOS)
        con reglas de mezclado de van der Waals

        Válida para altas temperaturas (hasta 2000 K) y presiones
        Más precisa que reglas lineales para gases densos
        """
        R = 8.314  # J/(mol·K)

        # Parámetros críticos y factores acéntricos
        params_criticos = {
            "CO2": {"Tc": 304.13, "Pc": 7.377e6, "omega": 0.2236, "MM": 0.04401},
            "H2O": {"Tc": 647.10, "Pc": 22.064e6, "omega": 0.3443, "MM": 0.01802},
            "CH4": {"Tc": 190.56, "Pc": 4.599e6, "omega": 0.0115, "MM": 0.01604},
            "O2": {"Tc": 154.58, "Pc": 5.043e6, "omega": 0.0222, "MM": 0.03200},
            "N2": {"Tc": 126.19, "Pc": 3.396e6, "omega": 0.0377, "MM": 0.02801},
            "CO": {"Tc": 132.86, "Pc": 3.494e6, "omega": 0.0497, "MM": 0.02801},
            "H2": {"Tc": 33.19, "Pc": 1.313e6, "omega": -0.2160, "MM": 0.00202},
            "C2H6": {"Tc": 305.32, "Pc": 4.872e6, "omega": 0.0995, "MM": 0.03007},
            "C3H8": {"Tc": 369.83, "Pc": 4.248e6, "omega": 0.1523, "MM": 0.04410},
            "C2H5OH": {"Tc": 513.92, "Pc": 6.148e6, "omega": 0.6450, "MM": 0.04607}
        }

        # Parámetros de interacción binaria kij (Peng-Robinson)
        # Fuente: Søreide & Whitson (1992), Li & Firoozabadi (2009)
        kij_matrix = {
            ("CO2", "H2O"): 0.1896,  # Par no-ideal crítico (alta polaridad diferencial)
            ("H2O", "CO2"): 0.1896,  # Simétrico
            # Otros pares asumen kij = 0 (comportamiento cuasi-ideal)
        }

        # Calcular parámetros a y b para cada componente
        a_i = {}
        b_i = {}
        MM_mezcla = 0

        for comp, frac in self.composicion.items():
            if frac > 1e-6 and comp in params_criticos:
                p = params_criticos[comp]
                Tc = p["Tc"]
                Pc = p["Pc"]
                omega = p["omega"]

                # Parámetro a de PR (dependiente de temperatura)
                Tr = self.T / Tc  # Temperatura reducida
                kappa = 0.37464 + 1.54226*omega - 0.26992*omega**2
                alpha = (1 + kappa*(1 - np.sqrt(Tr)))**2
                a = 0.45724 * (R**2 * Tc**2 / Pc) * alpha

                # Parámetro b de PR (constante)
                b = 0.07780 * (R * Tc / Pc)

                a_i[comp] = a
                b_i[comp] = b
                MM_mezcla += frac * p["MM"]

        # Reglas de mezclado de van der Waals
        a_mezcla = 0
        b_mezcla = 0

        for comp1, frac1 in self.composicion.items():
            if frac1 > 1e-6 and comp1 in a_i:
                b_mezcla += frac1 * b_i[comp1]

                for comp2, frac2 in self.composicion.items():
                    if frac2 > 1e-6 and comp2 in a_i:
                        # Regla de mezcla para a con parámetro de interacción binaria kij
                        kij = kij_matrix.get((comp1, comp2), 0.0)  # Default kij = 0 si no está definido
                        a_ij = (1 - kij) * np.sqrt(a_i[comp1] * a_i[comp2])
                        a_mezcla += frac1 * frac2 * a_ij

        # Resolver ecuación cúbica de PR para factor de compresibilidad Z
        A = a_mezcla * self.P / (R**2 * self.T**2)
        B = b_mezcla * self.P / (R * self.T)

        # Coeficientes: Z³ - (1-B)Z² + (A-3B²-2B)Z - (AB-B²-B³) = 0
        coef = [1, -(1 - B), A - 3*B**2 - 2*B, -(A*B - B**2 - B**3)]

        # Resolver cúbica y tomar raíz mayor (vapor)
        raices = np.roots(coef)
        raices_reales = raices[np.isreal(raices)].real
        Z = np.max(raices_reales) if len(raices_reales) > 0 else 1.0

        # Densidad molar
        rho_molar = self.P / (Z * R * self.T)
        self.rho = rho_molar * MM_mezcla  # kg/m³

        # Entalpía y entropía usando componentes puros + correcciones PR
        h_ideal = 0
        s_ideal = 0

        for comp, frac in self.composicion.items():
            if frac > 1e-6 and comp in params_criticos:
                try:
                    fluido = self._convertir_nombre_fluido(comp)
                    state = AbstractState("HEOS", fluido)
                    state.update(CP.PT_INPUTS, self.P, self.T)
                    h_ideal += frac * state.hmolar()
                    s_ideal += frac * state.smolar()
                except:
                    # Fallback a Cp constante
                    Cp_avg = 40000  # J/(mol·K)
                    h_ideal += frac * Cp_avg * (self.T - 298.15)
                    s_ideal += frac * Cp_avg * np.log(self.T / 298.15)

        # Departure functions de PR
        sqrt2 = np.sqrt(2)
        ln_term = np.log((Z + (1 + sqrt2)*B) / (Z + (1 - sqrt2)*B))

        h_dep = R * self.T * (Z - 1 - (a_mezcla / (2*sqrt2*b_mezcla*R*self.T)) * ln_term)
        s_dep = R * (np.log(Z - B) - (a_mezcla / (2*sqrt2*b_mezcla*R*self.T)) * ln_term)

        self.h = h_ideal + h_dep
        self.s = s_ideal + s_dep

        # Entropía de mezclado
        for comp, frac in self.composicion.items():
            if frac > 1e-6:
                self.s -= R * frac * np.log(frac)

    def _h_ideal_gas(self, T_ref=298.15):
        """Entalpía de gas ideal aproximada"""
        Cp = 35.0  # J/(mol·K) valor aproximado
        return Cp * (self.T - T_ref)

    def _s_ideal_gas(self, T_ref=298.15, P_ref=101325):
        """Entropía de gas ideal aproximada"""
        Cp = 35.0  # J/(mol·K)
        R = 8.314
        return Cp * np.log(self.T / T_ref) - R * np.log(self.P / P_ref)


class Combustible:
    """Clase para manejar diferentes tipos de combustibles"""
    def __init__(self, tipo):
        self.tipo = tipo
        self.composicion = self._get_composicion()
        self.LHV = self._get_LHV()  # Lower Heating Value en J/mol

    def _get_composicion(self):
        """Retorna la composición molar del combustible"""
        composiciones = {
            "Gas Natural": {"CH4": 0.9596, "C2H6": 0.0303, "C3H8": 0.0101},
            "Gas de Síntesis": {"CO": 0.40, "H2": 0.50, "CH4": 0.05, "CO2": 0.05},
            "Propano": {"C3H8": 1.0},
            "Etanol": {"C2H5OH": 1.0}
        }
        return composiciones.get(self.tipo, {"CH4": 1.0})

    def _get_LHV(self):
        """Retorna el poder calorífico inferior en J/mol"""
        LHVs = {
            "Gas Natural": 802000,  # J/mol aproximado
            "Gas de Síntesis": 250000,  # J/mol aproximado (mezcla CO+H2+CH4)
            "Propano": 2043000,
            "Etanol": 1277000
        }
        return LHVs.get(self.tipo, 802000)

    def calcular_O2_estequiometrico(self, n_combustible=1.0):
        """Calcula los moles de O2 necesarios para combustión estequiométrica"""
        O2_necesario = 0

        for componente, fraccion in self.composicion.items():
            n_comp = n_combustible * fraccion

            if componente == "CH4":
                O2_necesario += 2 * n_comp  # CH4 + 2O2 -> CO2 + 2H2O
            elif componente == "C2H6":
                O2_necesario += 3.5 * n_comp  # C2H6 + 3.5O2 -> 2CO2 + 3H2O
            elif componente == "C3H8":
                O2_necesario += 5 * n_comp  # C3H8 + 5O2 -> 3CO2 + 4H2O
            elif componente == "C2H5OH":
                O2_necesario += 3 * n_comp  # C2H5OH + 3O2 -> 2CO2 + 3H2O
            elif componente == "CO":
                O2_necesario += 0.5 * n_comp  # CO + 0.5O2 -> CO2
            elif componente == "H2":
                O2_necesario += 0.5 * n_comp  # H2 + 0.5O2 -> H2O
            # CO2 y H2O son inertes, no consumen O2

        return O2_necesario

    def calcular_productos_combustion(self, n_combustible=1.0):
        """Calcula la composición de productos de combustión"""
        n_CO2 = 0
        n_H2O = 0

        for componente, fraccion in self.composicion.items():
            n_comp = n_combustible * fraccion

            if componente == "CH4":
                n_CO2 += 1 * n_comp
                n_H2O += 2 * n_comp
            elif componente == "C2H6":
                n_CO2 += 2 * n_comp
                n_H2O += 3 * n_comp
            elif componente == "C3H8":
                n_CO2 += 3 * n_comp
                n_H2O += 4 * n_comp
            elif componente == "C2H5OH":
                n_CO2 += 2 * n_comp
                n_H2O += 3 * n_comp
            elif componente == "CO":
                n_CO2 += 1 * n_comp  # CO + 0.5O2 -> CO2
            elif componente == "H2":
                n_H2O += 1 * n_comp  # H2 + 0.5O2 -> H2O
            elif componente == "CO2":
                n_CO2 += 1 * n_comp  # CO2 pasa directo (inerte)
            elif componente == "H2O":
                n_H2O += 1 * n_comp  # H2O pasa directo (inerte)

        total = n_CO2 + n_H2O
        if total > 0:
            return {"CO2": n_CO2/total, "H2O": n_H2O/total}, n_CO2 + n_H2O
        return {"CO2": 0.5, "H2O": 0.5}, 0


class SimuladorBrayton:
    """Clase principal del simulador del ciclo de Brayton"""
    def __init__(self, parametros):
        self.params = parametros
        self.corrientes = {}
        self.W_neto = 0
        self.Q_combustion = 0
        self.eta_cycle = 0  # Eficiencia del ciclo (turbina - compresor CO2)
        self.eta_O2 = 0     # Eficiencia con ASU (incluye penalidad ASU)
        self.eta_CCS = 0    # Eficiencia global (incluye también compresor combustible)
        # NOTA: metodo_mezcla ya no es necesario como atributo
        # La clase Corriente selecciona automáticamente:
        # - HEOS para componentes puros
        # - Peng-Robinson con van der Waals para mezclas

    def simular(self):
        """Ejecuta la simulación completa del ciclo"""
        # Inicializar combustible
        combustible = Combustible(self.params['tipo_combustible'])

        # Flujo molar de referencia (base de cálculo)
        n_combustible = self.params['flujo_combustible']  # mol/s

        # ====================================================================
        # CORRIENTE 1: Combustible comprimido
        # ====================================================================
        # Calcular la temperatura de salida del compresor considerando eficiencia
        T_entrada_comp = self.params['T_ambiente']
        P_entrada_comp = 101325  # Pa (1 atm)
        P_salida_comp = self.params['P_combustion']
        eta_compresor_fuel = self.params['eta_compresor_fuel']

        # Crear corriente de entrada al compresor
        corriente_entrada_comp = Corriente("Combustible entrada compresor",
                                          T=T_entrada_comp, P=P_entrada_comp,
                                          composicion=combustible.composicion,
                                          flujo_molar=n_combustible)
        corriente_entrada_comp.calcular_propiedades()  # Auto: HEOS para puros, PR para mezclas

        # Temperatura de salida isentrópica (estimación usando relación politrópica)
        gamma = 1.3  # Aproximación para hidrocarburos
        T1_ideal = T_entrada_comp * (P_salida_comp / P_entrada_comp)**((gamma-1)/gamma)

        # Temperatura real considerando eficiencia
        T1_real = T_entrada_comp + (T1_ideal - T_entrada_comp) / eta_compresor_fuel

        self.corrientes[1] = Corriente(
            "Combustible comprimido",
            T=T1_real,
            P=self.params['P_combustion'],
            composicion=combustible.composicion,
            flujo_molar=n_combustible
        )
        self.corrientes[1].calcular_propiedades()  # Auto: HEOS para puros, PR para mezclas

        # ====================================================================
        # CORRIENTE 2: Oxígeno de ASU
        # ====================================================================
        n_O2 = combustible.calcular_O2_estequiometrico(n_combustible)
        self.corrientes[2] = Corriente(
            "Oxígeno de ASU",
            T=self.params['T_ambiente'] + 10,
            P=self.params['P_combustion'],
            composicion={"O2": 1.0},
            flujo_molar=n_O2
        )
        self.corrientes[2].calcular_propiedades()  # Auto: HEOS para puros, PR para mezclas

        # ====================================================================
        # ITERACIÓN PARA CONVERGENCIA DE CORRIENTE 11
        # ====================================================================
        # La corriente 11 (CO2 precalentado) tiene doble dependencia circular:
        # 1. Temperatura: C11 → C3 → C4 → recuperador → C11
        # 2. Flujo: n_CO2_recirc → C3 → ... → C9 → n_CO2_recirc
        # Resolvemos iterativamente ambas variables

        fraccion_recirculacion = self.params['fraccion_recirculacion']
        epsilon_rec = self.params['efectividad_recuperador']
        T_combustion = self.params['T_combustion']
        P_salida_turbina = self.params['P_salida_turbina']

        # ====================================================================
        # VALIDACIÓN DE FRACCIÓN DE RECIRCULACIÓN
        # ====================================================================
        if fraccion_recirculacion < 0.0:
            raise ValueError(
                f"❌ Error: Fracción de recirculación = {fraccion_recirculacion:.3f} < 0\n"
                "   La fracción debe ser no negativa (0 ≤ f < 1)"
            )

        if fraccion_recirculacion >= 1.0:
            raise ValueError(
                f"❌ Error: Fracción de recirculación = {fraccion_recirculacion:.3f} ≥ 1.0\n"
                "   Esto causa flujo infinito en el ciclo.\n"
                "   Matemáticamente: n_recirc = f/(1-f) × n_CO2_combustion → ∞ cuando f→1\n"
                "   La fracción debe ser estrictamente menor a 1 (f < 1.0)"
            )

        if fraccion_recirculacion > 0.97:
            print(f"⚠ ADVERTENCIA: Fracción de recirculación = {fraccion_recirculacion:.1%} es extremadamente alta (>97%)")
            print(f"   Flujo de recirculación = {fraccion_recirculacion/(1-fraccion_recirculacion):.1f} × n_CO2_combustion")
            print("   Esto puede causar:")
            print("   • Convergencia muy lenta (>15 iteraciones)")
            print("   • Inestabilidad numérica")
            print("   • Flujo másico muy alto en turbina y compresores")
            print("   Nota: Valores 95-97% son típicos en Ciclo Allam, pero requieren cuidado numérico.\n")

        elif fraccion_recirculacion > 0.90:
            print(f"ℹ Info: Fracción de recirculación = {fraccion_recirculacion:.1%} es alta (>90%)")
            print(f"   Flujo de recirculación = {fraccion_recirculacion/(1-fraccion_recirculacion):.1f} × n_CO2_combustion")
            print("   Esto es NORMAL en oxicombustión para controlar T_combustión.\n")

        # Parámetros de iteración
        MAX_ITER = 20
        TOL_T = 0.1  # K (tolerancia en temperatura)
        TOL_N = 0.001  # mol/s (tolerancia en flujo molar)

        # Estimación inicial de T11 y n_CO2_recirculado
        T10_inicial = self.params['T_ambiente'] + 15
        T4_inicial = T_combustion * (P_salida_turbina / self.params['P_combustion'])**((1.33-1)/1.33)
        T11_old = T10_inicial + epsilon_rec * (T4_inicial - T10_inicial)

        # Estimación inicial del flujo: basado en productos de combustión sin recirculación
        comp_productos_inicial, n_productos_combustion_inicial = combustible.calcular_productos_combustion(n_combustible)
        n_CO2_productos_inicial = n_productos_combustion_inicial * comp_productos_inicial.get("CO2", 0.333)
        n_CO2_recirculado_old = n_CO2_productos_inicial * fraccion_recirculacion

        # BUCLE DE ITERACIÓN PARA CONVERGENCIA
        # Ahora incluye convergencia de T_combustion calculada mediante balance energético
        T_combustion_old = T_combustion  # Estimación inicial

        for iter_count in range(MAX_ITER):
            # ====================================================================
            # CORRIENTE 3: Productos de combustión (entrada turbina)
            # ====================================================================
            # En el combustor se mezclan: Combustible (C1) + O2 (C2) + CO2 recirculado (C11)

            comp_productos, n_productos_combustion = combustible.calcular_productos_combustion(n_combustible)

            # Flujo total a la salida del combustor = productos + CO2 recirculado
            # Usamos el valor de iteración anterior (n_CO2_recirculado_old)
            n_CO2_productos = n_productos_combustion * comp_productos.get("CO2", 0.333)
            n_CO2_recirculado = n_CO2_recirculado_old  # Usar valor de iteración anterior

            # Flujo total de gases en corriente 3
            n_total_productos = n_productos_combustion + n_CO2_recirculado

            # Composición de la corriente 3 (productos + CO2 recirculado)
            if n_total_productos > 0:
                x_CO2_total = (n_CO2_productos + n_CO2_recirculado) / n_total_productos
                x_H2O_total = (n_productos_combustion * comp_productos.get("H2O", 0.667)) / n_total_productos
                comp_productos_total = {"CO2": x_CO2_total, "H2O": x_H2O_total}
            else:
                comp_productos_total = comp_productos

            # ====================================================================
            # CÁLCULO DE T_COMBUSTION MEDIANTE BALANCE ENERGÉTICO
            # ====================================================================
            # Balance CORREGIDO con estados de referencia:
            # H_productos(T) - H_productos(298K) = H_reactantes(T_in) - H_reactantes(298K) + Q_combustion
            #
            # Donde:
            # - H(T) son entalpías sensibles desde 298K (lo que calcula CoolProp)
            # - Q_combustion es el LHV (incluye diferencia de entalpía de formación)

            # Calcular DIFERENCIA de entalpía sensible entrada: ΔH = H(T_in) - H(298K)
            # Esto elimina dependencia de referencias arbitrarias de CoolProp

            delta_H_sensible_entrada = 0

            # Combustible: ΔH = n × [h(T_in) - h(298K)]
            h_fuel_Tin = self.corrientes[1].h if self.corrientes[1].h else 0
            corriente_ref_fuel_298 = Corriente("ref_fuel_298K", T=298.15, P=self.params['P_combustion'],
                                              composicion=combustible.composicion, flujo_molar=1.0)
            corriente_ref_fuel_298.calcular_propiedades()
            h_fuel_298 = corriente_ref_fuel_298.h if corriente_ref_fuel_298.h else 0
            delta_H_sensible_entrada += n_combustible * (h_fuel_Tin - h_fuel_298)

            # O2: ΔH = n × [h(T_in) - h(298K)]
            h_O2_Tin = self.corrientes[2].h if self.corrientes[2].h else 0
            corriente_ref_O2_298 = Corriente("ref_O2_298K", T=298.15, P=self.params['P_combustion'],
                                            composicion={"O2": 1.0}, flujo_molar=1.0)
            corriente_ref_O2_298.calcular_propiedades()
            h_O2_298 = corriente_ref_O2_298.h if corriente_ref_O2_298.h else 0
            delta_H_sensible_entrada += n_O2 * (h_O2_Tin - h_O2_298)

            # CO2 recirculado: ΔH = n × [h(T_in) - h(298K)]
            if 11 in self.corrientes and self.corrientes[11].h:
                h_CO2_Tin = self.corrientes[11].h
            else:
                # Estimación primera iteración: CO2 a ~968°C
                T_CO2_recirc_estimado = 1241  # K
                corriente_CO2_est = Corriente("CO2_recirc_estimado",
                                              T=T_CO2_recirc_estimado,
                                              P=self.params['P_combustion'],
                                              composicion={"CO2": 1.0}, flujo_molar=1.0)
                corriente_CO2_est.calcular_propiedades()
                h_CO2_Tin = corriente_CO2_est.h if corriente_CO2_est.h else 0

            corriente_ref_CO2_298 = Corriente("ref_CO2_298K", T=298.15, P=self.params['P_combustion'],
                                             composicion={"CO2": 1.0}, flujo_molar=1.0)
            corriente_ref_CO2_298.calcular_propiedades()
            h_CO2_298 = corriente_ref_CO2_298.h if corriente_ref_CO2_298.h else 0
            delta_H_sensible_entrada += n_CO2_recirculado * (h_CO2_Tin - h_CO2_298)

            # ============================================================================
            # BALANCE ENERGÉTICO SIMPLE (MÉTODO ORIGINAL RESTAURADO)
            # ============================================================================
            #
            # Balance: H_entrada + Q_combustion = H_salida
            #
            # Donde:
            # - H_entrada = suma de entalpías de C1, C2, C11 (J/s)
            # - Q_combustion = n_combustible × LHV (J/s)
            # - H_salida = n_productos × h_productos(T_combustión) (J/s)
            #
            # Calculamos H_entrada usando entalpías absolutas de CoolProp a P_combustion:

            H_entrada_combustor = 0

            # C1: Combustible
            if self.corrientes[1].h:
                H_entrada_combustor += n_combustible * self.corrientes[1].h

            # C2: Oxígeno
            if self.corrientes[2].h:
                H_entrada_combustor += n_O2 * self.corrientes[2].h

            # C11: CO2 recirculado (si existe)
            if 11 in self.corrientes and self.corrientes[11].h:
                H_entrada_combustor += n_CO2_recirculado * self.corrientes[11].h
            else:
                # Estimación primera iteración
                T_CO2_est = 1241  # K
                corr_est = Corriente("CO2_est", T=T_CO2_est, P=self.params['P_combustion'],
                                    composicion={"CO2": 1.0}, flujo_molar=n_CO2_recirculado)
                corr_est.calcular_propiedades()
                if corr_est.h:
                    H_entrada_combustor += n_CO2_recirculado * corr_est.h

            # ============================================================================
            # CORRECCIÓN DEL LHV PARA ALTA PRESIÓN
            # ============================================================================
            # El LHV está definido a 1 atm y 298K:
            #   LHV = H_productos(298K, 1atm) - H_reactivos(298K, 1atm)
            #
            # Para usar en el balance a P_combustion (150 bar), debemos corregir:
            #   LHV(P) = LHV(1atm) + ΔH_productos(P) - ΔH_reactivos(P)
            #
            # Donde:
            #   ΔH_productos(P) = H_productos(298K, P) - H_productos(298K, 1atm)
            #   ΔH_reactivos(P) = H_reactivos(298K, P) - H_reactivos(298K, 1atm)

            P_ref = 101325  # 1 atm
            T_ref = 298.15  # K
            P_combustion = self.params['P_combustion']

            # 1. Corrección para PRODUCTOS (CO2 + H2O)
            # Calcular productos de combustión sin recirculación
            comp_productos_puros, n_productos_combustion = combustible.calcular_productos_combustion(1.0)

            Delta_H_productos = 0
            for componente, fraccion in comp_productos_puros.items():
                if fraccion > 1e-6:
                    try:
                        # h(298K, 150bar)
                        corr_P = Corriente(f"{componente}_P", T=T_ref, P=P_combustion,
                                          composicion={componente: 1.0}, flujo_molar=1.0)
                        corr_P.calcular_propiedades()

                        # h(298K, 1atm)
                        corr_1atm = Corriente(f"{componente}_1atm", T=T_ref, P=P_ref,
                                             composicion={componente: 1.0}, flujo_molar=1.0)
                        corr_1atm.calcular_propiedades()

                        if corr_P.h and corr_1atm.h:
                            Delta_H_productos += fraccion * (corr_P.h - corr_1atm.h)
                    except:
                        pass  # Si falla, asumir corrección = 0 para ese componente

            # 2. Corrección para REACTIVOS (combustible + O2)
            Delta_H_reactivos = 0

            # Combustible
            for componente, fraccion in combustible.composicion.items():
                if fraccion > 1e-6:
                    try:
                        corr_P = Corriente(f"fuel_{componente}_P", T=T_ref, P=P_combustion,
                                          composicion={componente: 1.0}, flujo_molar=1.0)
                        corr_P.calcular_propiedades()

                        corr_1atm = Corriente(f"fuel_{componente}_1atm", T=T_ref, P=P_ref,
                                             composicion={componente: 1.0}, flujo_molar=1.0)
                        corr_1atm.calcular_propiedades()

                        if corr_P.h and corr_1atm.h:
                            Delta_H_reactivos += fraccion * (corr_P.h - corr_1atm.h)
                    except:
                        pass

            # O2
            n_O2_por_fuel = combustible.calcular_O2_estequiometrico(1.0)
            n_total_reactivos = 1.0 + n_O2_por_fuel
            frac_O2 = n_O2_por_fuel / n_total_reactivos

            try:
                corr_O2_P = Corriente("O2_P", T=T_ref, P=P_combustion,
                                     composicion={"O2": 1.0}, flujo_molar=1.0)
                corr_O2_P.calcular_propiedades()

                corr_O2_1atm = Corriente("O2_1atm", T=T_ref, P=P_ref,
                                        composicion={"O2": 1.0}, flujo_molar=1.0)
                corr_O2_1atm.calcular_propiedades()

                if corr_O2_P.h and corr_O2_1atm.h:
                    Delta_H_reactivos += frac_O2 * n_total_reactivos * (corr_O2_P.h - corr_O2_1atm.h)
            except:
                pass

            # LHV corregido para presión alta
            LHV_corregido = combustible.LHV + (Delta_H_productos - Delta_H_reactivos)

            # DEBUG: Mostrar correcciones
            if not hasattr(self, '_debug_LHV_printed'):
                print(f"\n{'='*80}")
                print(f"CORRECCIÓN DE LHV PARA ALTA PRESIÓN")
                print(f"{'='*80}")
                print(f"LHV(1 atm):              {combustible.LHV/1000:.2f} kJ/mol")
                print(f"ΔH_productos(P):         {Delta_H_productos/1000:+.2f} kJ/mol")
                print(f"ΔH_reactivos(P):         {Delta_H_reactivos/1000:+.2f} kJ/mol")
                print(f"Corrección neta:         {(Delta_H_productos - Delta_H_reactivos)/1000:+.2f} kJ/mol")
                print(f"LHV(150 bar):            {LHV_corregido/1000:.2f} kJ/mol")
                print(f"Cambio relativo:         {((LHV_corregido/combustible.LHV - 1)*100):+.2f}%")
                print(f"{'='*80}\n")
                self._debug_LHV_printed = True

            # Energía química con LHV corregido
            Q_combustion = n_combustible * LHV_corregido  # J/s

            # El balance se resuelve buscando T tal que:
            # H_productos(T) - H_productos(298K) = delta_H_sensible_entrada + LHV
            #
            # Es decir, la función objetivo pasa delta_H_sensible_entrada y LHV
            # y busca T tal que delta_H_sensible_productos = delta_H_sensible_entrada + LHV

            # DEBUG: Mostrar balance simple
            if not hasattr(self, '_debug_printed'):
                print(f"\n{'='*80}")
                print(f"COMBUSTIÓN - MÉTODO ORIGINAL SIMPLE (RESTAURADO)")
                print(f"{'='*80}")
                print(f"n_combustible:      {n_combustible:.2f} mol/s")
                print(f"n_O2:               {n_O2:.2f} mol/s")
                print(f"n_CO2_recirculado:  {n_CO2_recirculado:.2f} mol/s")
                print(f"n_total_productos:  {n_total_productos:.2f} mol/s")
                print(f"")
                print(f"H_entrada:          {H_entrada_combustor/1e6:.4f} MW")
                print(f"Q_combustion (LHV): {Q_combustion/1e6:.4f} MW")
                print(f"H_objetivo:         {(H_entrada_combustor + Q_combustion)/1e6:.4f} MW")
                print(f"{'='*80}\n")
                self._debug_printed = True

            # Resolver para T_combustion usando balance simple original:
            # H_salida = H_entrada + Q_combustion
            T_combustion_calculada = calcular_T_combustion_balance(
                H_entrada_combustor, n_total_productos, Q_combustion,
                self.params['P_combustion'], comp_productos_total
            )

            if hasattr(self, '_debug_printed') and self._debug_printed:
                print(f"\n--- RESULTADO FINAL ---")
                print(f"T_combustion_calculada: {T_combustion_calculada-273.15:.2f} °C")
                print(f"Composición productos: {comp_productos_total}")
                print(f"{'='*80}\n")
                self._debug_printed = False  # Reset para próxima simulación

            # Crear corriente C3 con la T_combustion calculada
            self.corrientes[3] = Corriente(
                "Productos combustión (entrada turbina)",
                T=T_combustion_calculada,
                P=self.params['P_combustion'],
                composicion=comp_productos_total,
                flujo_molar=n_total_productos
            )
            # Usar HEOS con mezcla por fugacidad para combustión (más preciso)
            self.corrientes[3].calcular_propiedades(metodo_mezcla="fugacidad")

            # Guardar flujos y LHV corregido para cálculos posteriores
            self._n_CO2_recirculado = n_CO2_recirculado
            self._n_productos_combustion = n_productos_combustion
            self._LHV_corregido = LHV_corregido  # LHV ajustado para P_combustion

            # ====================================================================
            # CORRIENTE 4: Salida de turbina
            # ====================================================================
            # Expansión isentrópica con eficiencia
            eta_turbina = self.params['eta_turbina']

            # Temperatura de salida ideal (isentrópica)
            T4_ideal = self.corrientes[3].T * (P_salida_turbina / self.params['P_combustion'])**((1.33-1)/1.33)

            # Temperatura real considerando eficiencia
            T4_real = self.corrientes[3].T - eta_turbina * (self.corrientes[3].T - T4_ideal)

            self.corrientes[4] = Corriente(
                "Salida turbina",
                T=T4_real,
                P=P_salida_turbina,
                composicion=comp_productos_total,
                flujo_molar=n_total_productos
            )
            self.corrientes[4].calcular_propiedades()  # Auto: PR (mezcla)

            # ====================================================================
            # CALCULAR C5-C9 PARA OBTENER n_CO2_recirculado_new
            # ====================================================================

            # C5: Salida recuperador (lado caliente)
            T5 = self.corrientes[4].T - epsilon_rec * (self.corrientes[4].T - (self.params['T_ambiente'] + 50))

            # C6: Entrada separador
            T6 = self.params['T_separador']

            # C7: CO2 puro (después de separar agua)
            n_CO2_puro = n_total_productos * comp_productos_total.get("CO2", 0.5)

            # C8: CO2 comprimido (el flujo no cambia)
            # C9: CO2 enfriado (el flujo no cambia)

            # CALCULAR n_CO2_recirculado_new basado en C9
            n_CO2_recirculado_new = fraccion_recirculacion * n_CO2_puro

            # Calcular T11_new basada en T4 real
            T11_new = T10_inicial + epsilon_rec * (T4_real - T10_inicial)

            # Verificar convergencia de TRES variables: T11, n_CO2_recirc, T_combustion
            error_T11 = abs(T11_new - T11_old)
            error_n = abs(n_CO2_recirculado_new - n_CO2_recirculado_old)
            error_Tcomb = abs(T_combustion_calculada - T_combustion_old)

            if error_T11 < TOL_T and error_n < TOL_N and error_Tcomb < TOL_T:
                # Convergió
                if iter_count > 0:  # Solo mostrar si hubo iteración
                    print(f"✓ Convergió en {iter_count+1} iteraciones (error_T11 = {error_T11:.3f} K, error_n = {error_n:.4f} mol/s, error_Tcomb = {error_Tcomb:.3f} K)")
                break

            # Actualizar para siguiente iteración
            T11_old = T11_new
            n_CO2_recirculado_old = n_CO2_recirculado_new
            T_combustion_old = T_combustion_calculada
        else:
            # No convergió en MAX_ITER iteraciones
            print(f"⚠ Advertencia: No convergió en {MAX_ITER} iteraciones (error_T11 = {error_T11:.3f} K, error_n = {error_n:.4f} mol/s, error_Tcomb = {error_Tcomb:.3f} K)")

        # Continuar con el resto de corrientes (C5-C12) fuera del bucle con valores convergidos

        # ====================================================================
        # CORRIENTE 5: Salida recuperador (enfriada)
        # ====================================================================
        # Efectividad del recuperador
        epsilon_rec = self.params['efectividad_recuperador']

        # Recirculación de CO2
        fraccion_recirculacion = self.params['fraccion_recirculacion']

        # Temperatura de salida del recuperador
        T5 = self.corrientes[4].T - epsilon_rec * (self.corrientes[4].T - (self.params['T_ambiente'] + 50))

        self.corrientes[5] = Corriente(
            "Salida recuperador",
            T=T5,
            P=P_salida_turbina * 0.98,  # Pequeña caída de presión
            composicion=comp_productos_total,
            flujo_molar=n_total_productos
        )
        self.corrientes[5].calcular_propiedades()  # Auto: PR (mezcla)

        # ====================================================================
        # CORRIENTE 6: Entrada separador agua
        # ====================================================================
        T6 = self.params['T_separador']

        self.corrientes[6] = Corriente(
            "Entrada separador",
            T=T6,
            P=P_salida_turbina * 0.96,
            composicion=comp_productos_total,
            flujo_molar=n_total_productos
        )
        self.corrientes[6].calcular_propiedades()  # Auto: PR (mezcla con H2O)

        # ====================================================================
        # CORRIENTE 7: CO2 puro (salida separador agua) - ENTRA A SEPARADOR DE FLUJOS
        # ====================================================================
        # El flujo de CO2 puro es todo el CO2 de los productos de combustión + recirculado
        n_CO2_puro_total = n_total_productos * comp_productos_total.get("CO2", 0.5)

        self.corrientes[7] = Corriente(
            "CO2 puro (salida separador agua)",
            T=T6,
            P=P_salida_turbina * 0.95,
            composicion={"CO2": 1.0},
            flujo_molar=n_CO2_puro_total
        )
        self.corrientes[7].calcular_propiedades()  # Auto: HEOS (CO2 puro)

        # ====================================================================
        # DIVISIÓN: CO2 capturado vs CO2 recirculado (DIAGRAMA CORREGIDO)
        # ====================================================================
        # La división ocurre DESPUÉS del separador de agua (C7)
        # Solo se comprime el CO2 que se recircula (C8), NO todo el CO2

        # Calcular n_CO2_recirculado basado en iteración
        n_CO2_recirculado_actual = n_CO2_recirculado_old  # De iteración anterior
        n_CO2_capturado = n_CO2_puro_total - n_CO2_recirculado_actual

        # ====================================================================
        # CORRIENTE 8: CO2 para recirculación (salida separador) → ENTRA AL COMPRESOR
        # ====================================================================
        self.corrientes[8] = Corriente(
            "CO2 para recirculación (entrada compresor)",
            T=T6,
            P=P_salida_turbina * 0.95,
            composicion={"CO2": 1.0},
            flujo_molar=n_CO2_recirculado_actual
        )
        self.corrientes[8].calcular_propiedades()  # Auto: HEOS (CO2 puro)

        # ====================================================================
        # CORRIENTE 12: CO2 capturado → VA A ALMACENAMIENTO (SIN COMPRIMIR)
        # ====================================================================
        # DIAGRAMA CORREGIDO: El CO2 capturado sale a P_salida_turbina, no se comprime
        self.corrientes[12] = Corriente(
            "CO2 capturado (a almacenamiento)",
            T=T6,
            P=P_salida_turbina * 0.95,  # Baja presión
            composicion={"CO2": 1.0},
            flujo_molar=n_CO2_capturado
        )
        self.corrientes[12].calcular_propiedades()  # Auto: HEOS (CO2 puro)

        # ====================================================================
        # CORRIENTE 9: CO2 comprimido → SALE DEL COMPRESOR, ENTRA A INTERCAMBIADOR 3
        # ====================================================================
        P_recirculacion = self.params['P_recirculacion']
        eta_compresor_CO2 = self.params['eta_compresor_CO2']

        # Compresión con eficiencia (C8 → C9)
        T9_ideal = self.corrientes[8].T * (P_recirculacion / self.corrientes[8].P)**((1.33-1)/1.33)
        T9_real = self.corrientes[8].T + (T9_ideal - self.corrientes[8].T) / eta_compresor_CO2

        self.corrientes[9] = Corriente(
            "CO2 comprimido (salida compresor)",
            T=T9_real,
            P=P_recirculacion,
            composicion={"CO2": 1.0},
            flujo_molar=n_CO2_recirculado_actual
        )
        self.corrientes[9].calcular_propiedades()  # Auto: HEOS (CO2 puro)

        # ====================================================================
        # CORRIENTE 10: CO2 enfriado → SALE INTERCAMBIADOR 3, ENTRA AL RECUPERADOR
        # ====================================================================
        T10 = self.params['T_ambiente'] + 15

        self.corrientes[10] = Corriente(
            "CO2 enfriado (salida intercambiador 3)",
            T=T10,
            P=P_recirculacion * 0.98,
            composicion={"CO2": 1.0},
            flujo_molar=n_CO2_recirculado_actual
        )
        self.corrientes[10].calcular_propiedades()  # Auto: HEOS (CO2 puro)

        # ====================================================================
        # CORRIENTE 11: CO2 precalentado → SALE RECUPERADOR, ENTRA A COMBUSTIÓN
        # ====================================================================
        # Calentamiento en el recuperador (C10 → C11)
        T11 = self.corrientes[10].T + epsilon_rec * (self.corrientes[4].T - self.corrientes[10].T)

        self.corrientes[11] = Corriente(
            "CO2 precalentado (salida recuperador)",
            T=T11,
            P=self.params['P_combustion'] * 0.95,
            composicion={"CO2": 1.0},
            flujo_molar=n_CO2_recirculado_actual
        )
        self.corrientes[11].calcular_propiedades()  # Auto: HEOS (CO2 puro)

        # ====================================================================
        # CÁLCULOS ENERGÉTICOS
        # ====================================================================
        # Siguiendo la notación de las fórmulas proporcionadas:
        # W_turb = m3 * (h3 - h4)
        # W_comp = m_combustible * (h1 - h_entrada)
        # eta_cycle = (W_turb - W_comp) / (mf * LHV)

        # Trabajo de la turbina: W_turb = m3 * (h3 - h4)
        # Donde m3 es el flujo de productos de combustión
        if self.corrientes[3].h and self.corrientes[4].h:
            W_turb = n_total_productos * (self.corrientes[3].h - self.corrientes[4].h) / 1e6  # MW
        else:
            W_turb = 0

        # Guardar como atributo para análisis de sensibilidad
        self.W_turbina = W_turb

        # Trabajo compresor combustible: W_comp = mf * (h1 - h_entrada)
        # Ya incluye el efecto de la eficiencia en la temperatura de salida
        if self.corrientes[1].h and corriente_entrada_comp.h:
            W_comp = n_combustible * (self.corrientes[1].h - corriente_entrada_comp.h) / 1e6  # MW
        else:
            # Aproximación si falla
            Cp_aprox = 35000  # J/(mol·K) aproximado
            Delta_T = T1_real - T_entrada_comp
            W_comp = n_combustible * Cp_aprox * Delta_T / 1e6  # MW

        # Trabajo compresor CO2: W_CO2comp = mCO2_recirc * (h9 - h8)
        # DIAGRAMA CORREGIDO: La división ocurre DESPUÉS del separador de agua (C7)
        # El Compresor 2 (C8→C9) actúa SOLO sobre la fracción recirculada (C8).
        # Por lo tanto, W_CO2comp es MUCHO MENOR que en la configuración anterior.
        if 8 in self.corrientes and 9 in self.corrientes and self.corrientes[8].h and self.corrientes[9].h:
            # C8 → C9: entrada y salida del compresor
            W_CO2comp_recirculacion = n_CO2_recirculado_actual * (self.corrientes[9].h - self.corrientes[8].h) / 1e6  # MW
            # En la nueva configuración, solo se comprime la fracción recirculada
            W_CO2comp_total = W_CO2comp_recirculacion  # Ya no hay "total", solo recirculación
        else:
            W_CO2comp_total = 0
            W_CO2comp_recirculacion = 0

        # Trabajo ASU: W_ASU = mO2 * e_ASU
        # Valor típico: 200-250 kWh/ton O2 = ~7000 J/mol O2
        W_ASU = n_O2 * 7000 / 1e6  # MW

        # ====================================================================
        # APORTE TÉRMICO EFECTIVO considerando el recuperador
        # ====================================================================
        # Balance energético del combustor (método 1: entalpía de entrada):
        # H_entrada = H_fuel + H_O2 + H_CO2_recirc
        # Q_necesario = H_salida - H_entrada
        #
        # El CO2 precalentado aporta entalpía, reduciendo el combustible necesario
        # para alcanzar T_combustion

        # Entalpía total de entrada al combustor
        H_entrada_combustor = 0

        # Contribución del combustible (C1)
        if self.corrientes[1].h:
            H_entrada_combustor += n_combustible * self.corrientes[1].h

        # Contribución del O2 (C2)
        if self.corrientes[2].h:
            H_entrada_combustor += n_O2 * self.corrientes[2].h

        # Contribución del CO2 recirculado (C11) - esta es la clave del recuperador
        if self.corrientes[11].h:
            H_entrada_combustor += n_CO2_recirculado * self.corrientes[11].h

        # Entalpía de salida del combustor (C3)
        if self.corrientes[3].h:
            H_salida_combustor = n_total_productos * self.corrientes[3].h
        else:
            H_salida_combustor = H_entrada_combustor

        # ====================================================================
        # CALOR DE ENTRADA PARA CÁLCULO DE EFICIENCIA
        # ====================================================================
        # IMPORTANTE: La eficiencia térmica debe calcularse con el calor REAL aportado
        # por el combustible, NO con el calor neto del combustor.
        #
        # Q_in = n_combustible × LHV (energía química liberada en la combustión)
        #
        # El recuperador NO aporta energía nueva, solo redistribuye la energía
        # ya generada mediante precalentamiento del CO2 recirculado.
        #
        # Usar LHV corregido a la presión del proceso si está disponible
        LHV_para_eficiencia = self._LHV_corregido if hasattr(self, '_LHV_corregido') else combustible.LHV
        Q_in_eficiencia = n_combustible * LHV_para_eficiencia / 1e6  # MW

        # Guardar para visualización (mismo valor que Q_in_eficiencia)
        self.Q_combustion = Q_in_eficiencia

        # BENEFICIO DEL RECUPERADOR:
        # - El CO2 precalentado (C11) reduce la necesidad de combustible para alcanzar T_combustion
        # - Esto se traduce en MAYOR trabajo neto (W_turb aumenta más que W_CO2comp)
        # - La eficiencia mejora porque W_neto aumenta, NO porque Q_in disminuye

        # Trabajo neto del sistema completo (considera TODO el trabajo de compresión de CO2)
        self.W_neto = W_turb - W_comp - W_CO2comp_recirculacion - W_ASU

        # Guardar trabajos individuales para visualización
        self.W_turbina = W_turb
        self.W_compresor_fuel = W_comp
        self.W_CO2comp_total = W_CO2comp_total
        self.W_CO2comp_recirculacion = W_CO2comp_recirculacion
        self.W_ASU = W_ASU

        # Guardar T_combustion calculada (después de convergencia)
        self.T_combustion_calculada = self.corrientes[3].T if 3 in self.corrientes else T_combustion

        # Eficiencias térmicas según las fórmulas proporcionadas
        if Q_in_eficiencia > 0:
            # 1. Eficiencia del ciclo: Turbina - Compresor CO2 recirculación
            # eta_cycle = (W_turb - W_comp_CO2_recirculacion) / (n_combustible × LHV)
            self.eta_cycle = ((W_turb - W_CO2comp_recirculacion) / Q_in_eficiencia) * 100

            # 2. Eficiencia con ASU: Incluye penalidad por producción de O2
            # eta_O2 = (W_turb - W_comp_CO2_recirculacion - W_ASU) / (n_combustible × LHV)
            self.eta_O2 = ((W_turb - W_CO2comp_recirculacion - W_ASU) / Q_in_eficiencia) * 100

            # 3. Eficiencia global: Incluye también compresor de combustible
            # eta_CCS = (W_turb - W_comp_CO2_recirculacion - W_ASU - W_comp_combustible) / (n_combustible × LHV)
            self.eta_CCS = ((W_turb - W_CO2comp_recirculacion - W_ASU - W_comp) / Q_in_eficiencia) * 100

            # DEBUG: Mostrar cálculo de eficiencias si hay valores anómalos
            if self.eta_cycle > 100 or self.eta_O2 > 100 or self.eta_CCS > 100:
                print(f"\n⚠️ ADVERTENCIA: Eficiencia >100% detectada")
                print(f"{'='*80}")
                print(f"ANÁLISIS DE EFICIENCIAS:")
                print(f"  W_turb:           {W_turb:.3f} MW")
                print(f"  W_comp_fuel:      {W_comp:.3f} MW")
                print(f"  W_ASU:            {W_ASU:.3f} MW")
                print(f"  W_CO2comp_recirc: {W_CO2comp_recirculacion:.3f} MW")
                print(f"  Q_in (LHV):       {Q_in_eficiencia:.3f} MW")
                print(f"  n_combustible:    {n_combustible:.2f} mol/s")
                print(f"  n_total_productos:{n_total_productos:.2f} mol/s")
                print(f"  Ratio n_prod/n_fuel: {n_total_productos/n_combustible:.1f}")
                print(f"")
                print(f"  W_neto = W_turb - W_CO2comp - W_ASU - W_comp_fuel")
                print(f"         = {W_turb:.3f} - {W_CO2comp_recirculacion:.3f} - {W_ASU:.3f} - {W_comp:.3f}")
                print(f"         = {W_turb - W_CO2comp_recirculacion - W_ASU - W_comp:.3f} MW")
                print(f"")
                print(f"  eta_cycle = (W_turb - W_CO2comp) / Q_in")
                print(f"            = ({W_turb:.3f} - {W_CO2comp_recirculacion:.3f}) / {Q_in_eficiencia:.3f}")
                print(f"            = {self.eta_cycle:.2f}%")
                print(f"  eta_O2 = (W_turb - W_CO2comp - W_ASU) / Q_in")
                print(f"         = ({W_turb:.3f} - {W_CO2comp_recirculacion:.3f} - {W_ASU:.3f}) / {Q_in_eficiencia:.3f}")
                print(f"         = {self.eta_O2:.2f}%")
                print(f"  eta_CCS = (W_turb - W_CO2comp - W_ASU - W_comp_fuel) / Q_in")
                print(f"          = ({W_turb:.3f} - {W_CO2comp_recirculacion:.3f} - {W_ASU:.3f} - {W_comp:.3f}) / {Q_in_eficiencia:.3f}")
                print(f"          = {self.eta_CCS:.2f}%")
                print(f"")
                if 3 in self.corrientes and 4 in self.corrientes:
                    print(f"  T_combustion: {self.corrientes[3].T-273.15:.1f}°C")
                    print(f"  T_salida_turb: {self.corrientes[4].T-273.15:.1f}°C")
                    print(f"  Delta_T_turb: {(self.corrientes[3].T - self.corrientes[4].T):.1f} K")
                    if self.corrientes[3].h and self.corrientes[4].h:
                        print(f"  h_entrada_turb: {self.corrientes[3].h/1e6:.4f} MJ/mol")
                        print(f"  h_salida_turb: {self.corrientes[4].h/1e6:.4f} MJ/mol")
                        print(f"  Delta_h_turb: {(self.corrientes[3].h - self.corrientes[4].h)/1e6:.4f} MJ/mol")
                print(f"{'='*80}\n")

            # No limitar artificialmente - dejar que muestre el valor real
            # Si hay errores, se verán en valores fuera de rango esperado
        else:
            self.eta_cycle = 0
            self.eta_O2 = 0
            self.eta_CCS = 0

        return True


# ============================================================================
# INTERFAZ DE STREAMLIT
# ============================================================================

# Botón de simulación (ubicado justo debajo del título)
simular_btn = st.button("🔥 SIMULAR", type="primary")

# ============================================================================
# SIDEBAR - Parámetros de entrada
# ============================================================================

st.sidebar.header("⚙️ Parámetros de Simulación")

# Tipo de combustible
tipo_combustible = st.sidebar.selectbox(
    "Tipo de Combustible",
    ["Gas Natural", "Gas de Síntesis", "Propano", "Etanol"]
)

# Método de cálculo de propiedades
# ============================================================================
# Método de cálculo (AUTOMÁTICO)
# ============================================================================
st.sidebar.info(
    "📐 **Ecuación de Estado**\n\n"
    "**Compuestos puros:**\n"
    "- HEOS\n\n"
    "**Mezclas:**\n"
    "- Peng-Robinson + van der Waals"
)

st.sidebar.subheader("Condiciones Operacionales")

# Temperatura ambiente (condición de frontera, NO optimizable)
T_ambiente = st.sidebar.number_input(
    "Temperatura Ambiente (°C)",
    min_value=0.0,
    max_value=50.0,
    value=25.0,
    step=1.0,
    help="Temperatura ambiente del sitio de instalación. Esta es una condición de frontera ambiental, no un parámetro de diseño optimizable."
)

# Presión de combustión (también define P_recirculacion por restricción termodinámica)
# IMPORTANTE: min_value = Pc_CO2 = 7.377 MPa
# Razón: Como P_recirculacion = P_combustion, el CO2 en C8-C11 debe estar supercrítico
# (NO porque C3 deba ser supercrítico, sino porque C8-C11 deben serlo para recirculación)
P_combustion = st.sidebar.number_input(
    "Presión de Combustión / Almacenamiento (MPa)",
    min_value=7.377,
    max_value=30.0,
    value=30.0,
    step=0.5,
    help="⚠️ Restricción: P ≥ 7.377 MPa (Pc del CO₂) para recirculación supercrítico. Mayor presión → Mayor ratio de expansión → Mayor eficiencia. Valor óptimo: 10-20 MPa"
)

# Valor inicial para la primera iteración (será recalculado automáticamente)
T_combustion = 1600.0  # °C (estimación inicial)

# Presión salida turbina
P_salida_turbina = st.sidebar.number_input(
    "Presión Salida Turbina (MPa)",
    min_value=0.1,
    max_value=5.0,
    value=0.12,
    step=0.01,
    help="Menor presión → Mayor expansión → Mayor trabajo. Mínimo: 0.1 MPa (1 bar, atmosférica). Máximo recomendado: 0.5 MPa para facilitar condensación de H₂O."
)

# Validaciones de P_salida_turbina
if P_salida_turbina >= P_combustion:
    st.sidebar.error("❌ Error: P_salida_turbina debe ser MENOR que P_combustion para que haya expansión en la turbina.")

# Verificar si T_separador resultante es razonable
from CoolProp.CoolProp import PropsSI
T_sat_H2O_at_P_salida = PropsSI('T', 'P', P_salida_turbina * 1e6, 'Q', 0, 'Water') - 273.15  # °C
T_separador_resultante = T_sat_H2O_at_P_salida - 10  # °C

if T_separador_resultante > 120:
    st.sidebar.warning(f"⚠️ P_salida_turbina = {P_salida_turbina:.2f} MPa resulta en T_separador ≈ {T_separador_resultante:.1f}°C (muy alta para condensación eficiente)")
elif T_separador_resultante < T_ambiente:
    st.sidebar.warning(f"⚠️ T_separador calculada ({T_separador_resultante:.1f}°C) es menor que T_ambiente ({T_ambiente:.1f}°C). Se requiere enfriamiento activo.")

st.sidebar.subheader("Eficiencias de Equipos")

# Eficiencia turbina
eta_turbina = st.sidebar.slider(
    "Eficiencia Isentrópica Turbina (%)",
    min_value=70.0,
    max_value=95.0,
    value=93.0,
    step=1.0,
    help="Turbinas modernas de gas: 88-92%. Valor óptimo: 90%"
) / 100

# Eficiencia compresor combustible
eta_compresor_fuel = st.sidebar.slider(
    "Eficiencia Isentrópica Compresor Combustible (%)",
    min_value=70.0,
    max_value=95.0,
    value=89.0,
    step=1.0,
    help="Compresor del combustible (C1). Compresores centrífugo/axial modernos: 85-90%. Valor óptimo: 88%"
) / 100

# Eficiencia compresor CO2
eta_compresor_CO2 = st.sidebar.slider(
    "Eficiencia Isentrópica Compresor CO₂ (%)",
    min_value=70.0,
    max_value=95.0,
    value=89.0,
    step=1.0,
    help="Compresor de CO₂ para recirculación (C7→C8). Compresión de CO₂: 82-88%. Valor típico: 85%"
) / 100

# Efectividad recuperador
efectividad_recuperador = st.sidebar.slider(
    "Efectividad Recuperador (%)",
    min_value=50.0,
    max_value=95.0,
    value=90.0,
    step=5.0,
    help="Recuperadores avanzados: 85-95%. Mayor recuperación → Menor consumo combustible"
) / 100

st.sidebar.subheader("Parámetros del Ciclo")

# Flujo de combustible
flujo_combustible = st.sidebar.number_input(
    "Flujo Molar Combustible (mol/s)",
    min_value=1.0,
    max_value=1000.0,
    value=12.0,  # Ajustado para T_combustion < 1800°C con alta recirculación (93%)
    step=1.0,
    help="⚙️ Flujo de combustible al combustor.\n\n"
         "📊 Rango típico: 10-15 mol/s (para T_combustión < 1800°C con alta recirculación)\n\n"
         "🎯 Objetivo: Mantener T_combustión entre 1200-1800°C\n\n"
         "⚠️ Si T_combustión > 1800°C → REDUCIR flujo o AUMENTAR recirculación\n"
         "⚠️ Si T_combustión < 1100°C → AUMENTAR flujo"
)

# Fracción de recirculación
fraccion_recirculacion = st.sidebar.slider(
    "Fracción de Recirculación CO2 (%)",
    min_value=0.0,
    max_value=99.0,
    value=93.0,
    step=1.0,
    help="Rango típico Ciclo Allam: 85-95%. Valores >97% pueden causar inestabilidad numérica."
) / 100

# Advertencia para valores muy altos de recirculación
if fraccion_recirculacion > 0.97:
    st.sidebar.error(
        f"🚨 Recirculación extremadamente alta ({fraccion_recirculacion*100:.0f}%): "
        "Puede causar inestabilidad numérica y no es práctica (captura de CO₂ <3%)."
    )
elif fraccion_recirculacion > 0.95:
    st.sidebar.warning(
        f"⚠️ Recirculación muy alta ({fraccion_recirculacion*100:.0f}%): "
        "La convergencia puede ser lenta y la captura de CO₂ es muy baja (<5%)."
    )

# NOTA: P_recirculacion = P_combustion (restricción termodinámica)
# El CO2 recirculado no puede cambiar de presión en el recuperador (solo intercambiador de calor)
# Por lo tanto, el compresor de CO2 debe comprimir directamente a P_combustion
P_recirculacion = P_combustion  # Igualdad termodinámica obligatoria

# T_separador se calcula automáticamente
T_separador_K = calcular_T_separador_automatica(P_salida_turbina * 1e6)  # En K
T_separador = T_separador_K - 273.15  # En °C

# ============================================================================
# Recopilar parámetros
# ============================================================================
parametros = {
    'tipo_combustible': tipo_combustible,
    # 'metodo_mezcla' eliminado: ahora es automático (HEOS para puros, PR para mezclas)
    'T_ambiente': T_ambiente + 273.15,  # Convertir a K
    'P_combustion': P_combustion * 1e6,  # Convertir de MPa a Pa
    'T_combustion': T_combustion + 273.15,  # Convertir a K
    'P_salida_turbina': P_salida_turbina * 1e6,  # Convertir de MPa a Pa
    'eta_turbina': eta_turbina,
    'eta_compresor_fuel': eta_compresor_fuel,
    'eta_compresor_CO2': eta_compresor_CO2,
    'efectividad_recuperador': efectividad_recuperador,
    'flujo_combustible': flujo_combustible,
    'fraccion_recirculacion': fraccion_recirculacion,
    'P_recirculacion': P_recirculacion * 1e6,  # Convertir de MPa a Pa (igual a P_combustion)
    'T_separador': T_separador_K  # Ya está en K
}

# ============================================================================
# Ejecutar simulación y guardar en session_state
# ============================================================================
# Inicializar session_state para la simulación si no existe
if 'simulador' not in st.session_state:
    st.session_state['simulador'] = None
if 'simulacion_exitosa' not in st.session_state:
    st.session_state['simulacion_exitosa'] = False

if simular_btn:
    with st.spinner('Ejecutando simulación...'):
        try:
            st.session_state['simulador'] = SimuladorBrayton(parametros)
            st.session_state['simulacion_exitosa'] = st.session_state['simulador'].simular()
            if st.session_state['simulacion_exitosa']:
                st.success('✅ Simulación completada exitosamente!')
            else:
                st.error('❌ La simulación no convergió correctamente')
        except Exception as e:
            st.error(f'❌ Error en la simulación: {str(e)}')
            st.session_state['simulacion_exitosa'] = False
            st.session_state['simulador'] = None

# Recuperar simulación desde session_state
simulador = st.session_state['simulador']
simulacion_exitosa = st.session_state['simulacion_exitosa']

# ============================================================================
# MENÚ PRINCIPAL - Navegación jerárquica de 2 niveles
# ============================================================================

# Definir estructura del menú
main_menu_names = [
    "📊 Datos Simulación",
    "🔬 Análisis de Sensibilidad",
    "🎯 Optimización"
]

# Submenús para "Datos Simulación"
sub_menu_datos = [
    "📊 Diagrama del Proceso",
    "📋 Tabla de Propiedades",
    "📈 Diagrama P-H",
    "📉 Diagrama T-S"
]

# Inicializar índices en session_state si no existen
if 'main_menu_index' not in st.session_state:
    st.session_state['main_menu_index'] = 0
if 'sub_menu_index' not in st.session_state:
    st.session_state['sub_menu_index'] = 0

# Menú principal usando directamente el índice
selected_main_index = st.radio(
    "Seleccionar sección:",
    options=range(len(main_menu_names)),
    index=st.session_state['main_menu_index'],
    format_func=lambda x: main_menu_names[x],
    horizontal=True,
    key='main_menu_selector',
    label_visibility='collapsed'
)

# Actualizar índice del menú principal
st.session_state['main_menu_index'] = selected_main_index
selected_main = main_menu_names[selected_main_index]

# Submenú para "Datos Simulación"
if selected_main_index == 0:  # "📊 Datos Simulación"
    selected_sub_index = st.radio(
        "Vista de datos:",
        options=range(len(sub_menu_datos)),
        index=st.session_state['sub_menu_index'],
        format_func=lambda x: sub_menu_datos[x],
        horizontal=True,
        key='sub_menu_selector',
        label_visibility='collapsed'
    )
    st.session_state['sub_menu_index'] = selected_sub_index
    selected_sub = sub_menu_datos[selected_sub_index]
else:
    selected_sub = sub_menu_datos[0]

st.markdown("---")  # Separador visual

# Crear variables de pestaña para compatibilidad con código existente
class TabContext:
    def __init__(self, is_active):
        self.is_active = is_active
    def __enter__(self):
        return self
    def __exit__(self, *args):
        pass

# Asignar pestañas según el menú activo
tab1 = TabContext(selected_main_index == 0 and selected_sub == "📊 Diagrama del Proceso")
tab2 = TabContext(selected_main_index == 0 and selected_sub == "📋 Tabla de Propiedades")
tab3 = TabContext(selected_main_index == 0 and selected_sub == "📈 Diagrama P-H")
tab4 = TabContext(selected_main_index == 0 and selected_sub == "📉 Diagrama T-S")
tab5 = TabContext(selected_main_index == 1)  # Análisis de Sensibilidad
tab6 = TabContext(selected_main_index == 2)  # Optimización

# ============================================================================
# TAB 1: Diagrama del Proceso
# ============================================================================
if tab1.is_active:
    st.header("Diagrama del Ciclo de Brayton con Oxicombustión")

    # Intentar cargar el diagrama
    import os
    import streamlit.components.v1 as components

    try:
        # Obtener la ruta absoluta del directorio del script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        diagrama_path_svg = os.path.join(script_dir, "diagramas brayton.svg")
        diagrama_path_png = os.path.join(script_dir, "diagrama_brayton.png")
        # Preferir el diagrama corregido si existe
        diagrama_path_png_corr = os.path.join(script_dir, "diagrama_brayton_corregido.png")
        if os.path.exists(diagrama_path_png_corr):
            diagrama_path_png = diagrama_path_png_corr

        # Cargar imagen PNG
        with st.expander("Ver Diagrama del Proceso", expanded=True):
            if os.path.exists(diagrama_path_png):
                imagen = Image.open(diagrama_path_png)
                st.image(imagen, caption="Diagrama del proceso", use_container_width=True)
            else:
                # Intentar con ruta relativa
                try:
                    imagen = Image.open("diagrama_brayton.png")
                    st.image(imagen, caption="Diagrama del proceso", use_container_width=True)
                except:
                    st.warning("⚠️ No se pudo cargar el diagrama.")
    except FileNotFoundError:
        st.warning("⚠️ No se pudo cargar el diagrama. Asegúrate de que 'diagramas brayton.svg' o 'diagrama_brayton.png' estén en el mismo directorio que el script.")
        st.info(f"Directorio actual: {os.getcwd()}")
    except Exception as e:
        st.error(f"Error al cargar el diagrama: {str(e)}")

    if simulacion_exitosa and simulador:
        st.subheader("Resultados Energéticos")

        # Primera fila: Trabajo, Calor y Temperatura de Combustión
        col1, col2, col3 = st.columns(3)

        with col1:
            st.metric("Trabajo Neto", f"{simulador.W_neto:.2f} MW")

        with col2:
            st.metric("Calor de Combustión", f"{simulador.Q_combustion:.2f} MW")

        with col3:
            T_comb_calc_C = simulador.T_combustion_calculada - 273.15 if hasattr(simulador, 'T_combustion_calculada') else 0
            st.metric("T Combustión (calculada)", f"{T_comb_calc_C:.1f} °C",
                     help="Calculada mediante balance energético del combustor")

        # Advertencias dinámicas basadas en T_combustion calculada
        if hasattr(simulador, 'T_combustion_calculada'):
            T_comb_K = simulador.T_combustion_calculada
            T_comb_C = T_comb_K - 273.15

            # Nivel crítico: fuera de límites absolutos
            if T_comb_K >= 2000 + 273.15 - 10:  # > 1990°C
                flujo_sugerido = flujo_combustible * 0.65
                st.error(
                    f"🚨 **CRÍTICO**: T_combustión = {T_comb_C:.1f}°C está en el límite máximo!\n\n"
                    f"**Acción requerida**: Reducir flujo de combustible a **~{flujo_sugerido:.0f} mol/s**\n\n"
                    f"⚠️ Resultados actuales NO son realistas (saturación en 2000°C)"
                )
            elif T_comb_K <= 1000 + 273.15 + 10:  # < 1010°C
                flujo_sugerido = flujo_combustible * 1.4
                st.error(
                    f"🚨 **CRÍTICO**: T_combustión = {T_comb_C:.1f}°C está en el límite mínimo!\n\n"
                    f"**Acción requerida**: Aumentar flujo de combustible a **~{flujo_sugerido:.0f} mol/s** o reducir recirculación\n\n"
                    f"⚠️ Combustión inestable a esta temperatura"
                )

            # Nivel advertencia: fuera del rango ideal pero dentro de límites
            elif T_comb_C > 1800:  # 1800-1990°C
                flujo_sugerido = flujo_combustible * 0.85
                st.warning(
                    f"⚠️ **Advertencia**: T_combustión = {T_comb_C:.1f}°C es alta (límite materiales ~1800°C)\n\n"
                    f"**Recomendación**: Reducir flujo a **~{flujo_sugerido:.0f} mol/s** para mayor realismo"
                )
            elif T_comb_C < 1200:  # 1010-1200°C
                flujo_sugerido = flujo_combustible * 1.15
                st.warning(
                    f"⚠️ **Advertencia**: T_combustión = {T_comb_C:.1f}°C es baja para operación eficiente\n\n"
                    f"**Recomendación**: Aumentar flujo a **~{flujo_sugerido:.0f} mol/s** para mejor eficiencia"
                )

            # Nivel info: rango óptimo
            elif 1400 <= T_comb_C <= 1700:
                st.success(f"✅ T_combustión = {T_comb_C:.1f}°C está en el rango óptimo (1400-1700°C)")
            elif 1200 <= T_comb_C < 1400 or 1700 < T_comb_C <= 1800:
                st.info(f"ℹ️ T_combustión = {T_comb_C:.1f}°C está en rango aceptable (1200-1800°C)")

        # Segunda fila: Tres eficiencias
        st.markdown("---")
        st.subheader("Eficiencias Térmicas")

        col3, col4, col5 = st.columns(3)

        with col3:
            st.metric("η Ciclo", f"{simulador.eta_cycle:.2f} %",
                     help="Eficiencia del ciclo: turbina - compresor CO₂ recirculación")

        with col4:
            st.metric("η con ASU", f"{simulador.eta_O2:.2f} %",
                     help="Eficiencia con ASU: incluye penalidad de la unidad de separación de aire")

        with col5:
            st.metric("η Global", f"{simulador.eta_CCS:.2f} %",
                     help="Eficiencia global: incluye también compresor de combustible")

        # Panel de diagnóstico detallado
        with st.expander("🔍 Diagnóstico Detallado de Cálculos"):
            st.markdown("### 📊 FLUJOS MOLARES")

            # Recalcular flujos para verificación
            n_combustible = simulador.params['flujo_combustible']
            combustible = Combustible(simulador.params['tipo_combustible'])
            n_O2 = combustible.calcular_O2_estequiometrico(n_combustible)
            _, n_productos_combustion = combustible.calcular_productos_combustion(n_combustible)

            # Flujos en corriente 3 (con recirculación)
            n_total_C3 = simulador.corrientes[3].flujo_molar if 3 in simulador.corrientes else 0
            n_CO2_recirculado = simulador._n_CO2_recirculado if hasattr(simulador, '_n_CO2_recirculado') else 0

            # Flujo de CO2 puro (C7)
            if 7 in simulador.corrientes:
                n_CO2_total = simulador.corrientes[7].flujo_molar
            else:
                n_CO2_total = n_total_C3 * simulador.corrientes[3].composicion.get("CO2", 0.5)

            st.write(f"**n_combustible** = {n_combustible:.2f} mol/s")
            st.write(f"**n_O2** = {n_O2:.2f} mol/s (estequiométrico)")
            st.write(f"**n_productos_combustión** = {n_productos_combustion:.2f} mol/s (CO₂ + H₂O de combustión)")
            st.write(f"**n_CO2_recirculado** = {n_CO2_recirculado:.2f} mol/s (f_recirc = {simulador.params['fraccion_recirculacion']*100:.1f}%)")
            st.write(f"**n_total_C3** (turbina) = {n_total_C3:.2f} mol/s = productos + CO₂ recirculado")
            st.write(f"**n_CO2_total** (C7, separado) = {n_CO2_total:.2f} mol/s")

            st.markdown("---")
            st.markdown("### ⚡ TRABAJOS (MW)")

            # Trabajo turbina
            if simulador.corrientes[3].h and simulador.corrientes[4].h:
                W_turb_calc = n_total_C3 * (simulador.corrientes[3].h - simulador.corrientes[4].h) / 1e6
                st.write(f"**W_turbina** = {n_total_C3:.2f} × ({simulador.corrientes[3].h/1000:.2f} - {simulador.corrientes[4].h/1000:.2f}) = **{W_turb_calc:.3f} MW**")

            # Trabajo compresor combustible
            if 1 in simulador.corrientes:
                corriente_entrada = Corriente("temp", T=simulador.params['T_ambiente'], P=101325,
                                             composicion=combustible.composicion, flujo_molar=n_combustible)
                corriente_entrada.calcular_propiedades()  # Auto
                if corriente_entrada.h and simulador.corrientes[1].h:
                    W_comp_calc = n_combustible * (simulador.corrientes[1].h - corriente_entrada.h) / 1e6
                    st.write(f"**W_comp_fuel** = {n_combustible:.2f} × ({simulador.corrientes[1].h/1000:.2f} - {corriente_entrada.h/1000:.2f}) = **{W_comp_calc:.3f} MW**")

            # Trabajo compresor CO2 - DIAGRAMA CORREGIDO
            # La división ocurre DESPUÉS del separador de agua (C7)
            # C7 → separador → C8 (recirculación) → Compresor → C9 (comprimido)
            #                → C12 (captura, sin comprimir)
            if 8 in simulador.corrientes and 9 in simulador.corrientes:
                if simulador.corrientes[8].h and simulador.corrientes[9].h:
                    # SOLO se comprime la fracción recirculada (C8 → C9)
                    W_CO2_recirc_calc = n_CO2_recirculado * (simulador.corrientes[9].h - simulador.corrientes[8].h) / 1e6
                    st.write(f"**W_comp_CO2_recirculación** = {n_CO2_recirculado:.2f} × ({simulador.corrientes[9].h/1000:.2f} - {simulador.corrientes[8].h/1000:.2f}) = **{W_CO2_recirc_calc:.3f} MW**")

                    # MOSTRAR RATIO DE COMPRESIÓN (C8 → C9)
                    P8 = simulador.corrientes[8].P / 1e6  # MPa
                    P9 = simulador.corrientes[9].P / 1e6  # MPa
                    RC = P9 / P8
                    st.warning(f"⚠️ **Ratio de compresión CO₂:** {RC:.1f} (de {P8:.3f} MPa a {P9:.2f} MPa)")

            # Trabajo ASU
            W_ASU_calc = n_O2 * 7000 / 1e6
            st.write(f"**W_ASU** = {n_O2:.2f} × 7.0 kJ/mol = **{W_ASU_calc:.3f} MW**")

            st.markdown("---")
            st.markdown("### 🔥 CALOR Y EFICIENCIAS")

            # Calor de entrada
            Q_in_calc = n_combustible * combustible.LHV / 1e6
            st.write(f"**Q_in** = {n_combustible:.2f} × {combustible.LHV/1000:.2f} kJ/mol = **{Q_in_calc:.3f} MW**")

            # Trabajo neto y eficiencias
            W_neto_calc = simulador.W_neto
            st.write(f"**W_neto** = W_turb - W_comp_fuel - W_ASU - W_comp_CO2_recirc = **{W_neto_calc:.3f} MW**")

            st.write(f"**η_cycle** = (W_turb - W_comp_CO2_recirc) / Q_in = **{simulador.eta_cycle:.2f}%**")
            st.write(f"**η_O2** = (W_turb - W_comp_CO2_recirc - W_ASU) / Q_in = **{simulador.eta_O2:.2f}%**")
            st.write(f"**η_CCS (Global)** = (W_turb - W_comp_CO2_recirc - W_ASU - W_comp_fuel) / Q_in = **{simulador.eta_CCS:.2f}%**")
            st.info("ℹ️ **Nota:** η_cycle considera turbina y compresor CO₂. η_O2 añade penalidad ASU. η_CCS (global) incluye también compresor de combustible.")

            st.markdown("---")
            st.markdown("### 🌡️ TEMPERATURAS Y PRESIONES CRÍTICAS")
            st.write(f"**C3 (entrada turbina):** T = {simulador.corrientes[3].T-273.15:.1f}°C, P = {simulador.corrientes[3].P/1e6:.2f} MPa")
            st.write(f"**C4 (salida turbina):** T = {simulador.corrientes[4].T-273.15:.1f}°C, P = {simulador.corrientes[4].P/1e6:.3f} MPa")
            st.write(f"**C7 (CO₂ puro, salida separador agua):** T = {simulador.corrientes[7].T-273.15:.1f}°C, P = {simulador.corrientes[7].P/1e6:.3f} MPa, Flujo = {simulador.corrientes[7].flujo_molar:.2f} mol/s")
            if 8 in simulador.corrientes:
                st.write(f"**C8 (CO₂ recirculación, entrada compresor):** T = {simulador.corrientes[8].T-273.15:.1f}°C, P = {simulador.corrientes[8].P/1e6:.3f} MPa, Flujo = {simulador.corrientes[8].flujo_molar:.2f} mol/s")
            if 9 in simulador.corrientes:
                st.write(f"**C9 (CO₂ comprimido, salida compresor):** T = {simulador.corrientes[9].T-273.15:.1f}°C, P = {simulador.corrientes[9].P/1e6:.2f} MPa")
            if 10 in simulador.corrientes:
                st.write(f"**C10 (CO₂ enfriado, salida intercambiador):** T = {simulador.corrientes[10].T-273.15:.1f}°C, P = {simulador.corrientes[10].P/1e6:.3f} MPa")
            if 11 in simulador.corrientes:
                st.write(f"**C11 (CO₂ precalentado, entrada combustión):** T = {simulador.corrientes[11].T-273.15:.1f}°C, P = {simulador.corrientes[11].P/1e6:.3f} MPa")
            if 12 in simulador.corrientes:
                st.write(f"**C12 (CO₂ capturado, a almacenamiento):** T = {simulador.corrientes[12].T-273.15:.1f}°C, P = {simulador.corrientes[12].P/1e6:.3f} MPa, Flujo = {simulador.corrientes[12].flujo_molar:.2f} mol/s")

# ============================================================================
# TAB 2: Tabla de Propiedades
# ============================================================================
if tab2.is_active:
    st.header("Propiedades Termodinámicas de las Corrientes")

    if simulacion_exitosa and simulador:
        # Crear tabla de propiedades
        data = []
        for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]:
            if i in simulador.corrientes:
                c = simulador.corrientes[i]
                data.append({
                    "Corriente": f"{i}",
                    "Nombre": c.nombre,
                    "T (°C)": f"{c.T - 273.15:.1f}" if c.T else "-",
                    "P (MPa)": f"{c.P / 1e6:.3f}" if c.P else "-",
                    "H (kJ/mol)": f"{c.h / 1000:.2f}" if c.h else "-",
                    "S (kJ/mol·K)": f"{c.s / 1000:.4f}" if c.s else "-",
                    "Flujo (mol/s)": f"{c.flujo_molar:.2f}" if c.flujo_molar else "-",
                    "Composición": ", ".join([f"{comp}: {frac:.3f}" for comp, frac in c.composicion.items()])
                })

        df = pd.DataFrame(data)
        st.dataframe(df, width="stretch", hide_index=True)

        # Opción de descarga
        csv = df.to_csv(index=False)
        st.download_button(
            label="📥 Descargar tabla como CSV",
            data=csv,
            file_name="propiedades_corrientes.csv",
            mime="text/csv"
        )
    else:
        st.info("Presiona el botón 'SIMULAR' para ver las propiedades de las corrientes.")

# ============================================================================
# TAB 3: Diagrama P-H
# ============================================================================
if tab3.is_active:
    st.header("Diagrama Presión-Entalpía")

    if simulacion_exitosa and simulador:
        # Corrientes del ciclo: 3 -> 4 -> 5 -> [separador] -> 7 -> 8 -> 9 -> 10 -> 11
        # NOTA: Eliminamos C6 porque el separador NO es un proceso termodinámico continuo
        corrientes_ciclo = [3, 4, 5, 7, 8, 9, 10, 11]

        H_values = []
        P_values = []
        labels = []
        nombres = []

        # Colores para el ciclo principal
        color_map = {
            3: '#d62728',   # Rojo - Productos combustión (entrada turbina)
            4: '#ff7f0e',   # Naranja - Salida turbina
            5: '#bcbd22',   # Amarillo-verde - Salida recuperador (mezcla CO2+H2O)
            7: '#9467bd',   # Púrpura - CO2 puro (después del separador)
            8: '#8c564b',   # Café - CO2 comprimido
            9: '#e377c2',   # Rosa - CO2 enfriado
            10: '#7f7f7f',  # Gris - CO2 recirculado
            11: '#2ca02c'   # Verde - CO2 precalentado
        }

        for i in corrientes_ciclo:
            if i in simulador.corrientes:
                c = simulador.corrientes[i]
                if c.h and c.P:
                    H_values.append(c.h / 1000)  # kJ/mol
                    P_values.append(c.P / 1e6)  # MPa
                    labels.append(f"{i}")
                    nombres.append(c.nombre)

        # Crear gráfico con Plotly
        fig = go.Figure()

        # Línea conectando C3->C4->C5 (ciclo de potencia con mezcla)
        if len(H_values) >= 3:
            fig.add_trace(go.Scatter(
                x=H_values[0:3],  # C3, C4, C5
                y=P_values[0:3],
                mode='lines+markers',
                line=dict(color='rgba(100, 100, 100, 0.5)', width=2, dash='solid'),
                marker=dict(size=14, line=dict(width=2, color='black')),
                showlegend=False,
                hoverinfo='skip'
            ))

        # Discontinuidad C5->C7 (SEPARADOR - cambio de composición)
        if len(H_values) >= 4:
            idx_c5 = 2  # C5 es el tercer elemento (índice 2)
            idx_c7 = 3  # C7 es el cuarto elemento (índice 3)
            fig.add_trace(go.Scatter(
                x=[H_values[idx_c5], H_values[idx_c7]],
                y=[P_values[idx_c5], P_values[idx_c7]],
                mode='lines',
                line=dict(color='rgba(255, 0, 0, 0.4)', width=2, dash='dot'),  # Rojo punteado
                showlegend=False,
                hoverinfo='skip'
            ))

        # Línea conectando C7->C8->C9->C10->C11 (ciclo de CO2 puro)
        if len(H_values) >= 4:
            fig.add_trace(go.Scatter(
                x=H_values[3:],  # C7, C8, C9, C10, C11
                y=P_values[3:],
                mode='lines+markers',
                line=dict(color='rgba(100, 100, 100, 0.5)', width=2, dash='solid'),
                marker=dict(size=14, line=dict(width=2, color='black')),
                showlegend=False,
                hoverinfo='skip'
            ))

        # Línea cerrando el ciclo (11->3)
        if len(H_values) >= 2:
            fig.add_trace(go.Scatter(
                x=[H_values[-1], H_values[0]],  # Del último (C11) al primero (C3)
                y=[P_values[-1], P_values[0]],
                mode='lines',
                line=dict(color='rgba(0, 128, 0, 0.6)', width=3, dash='dash'),  # Verde punteado
                showlegend=False,
                hoverinfo='skip'
            ))

        # Puntos individuales con colores y etiquetas
        for i, (h, p, label, nombre) in enumerate(zip(H_values, P_values, labels, nombres)):
            corriente_id = corrientes_ciclo[i]
            fig.add_trace(go.Scatter(
                x=[h],
                y=[p],
                mode='markers+text',
                marker=dict(size=16, color=color_map[corriente_id], line=dict(width=2, color='black')),
                text=f"<b>{label}</b>",
                textposition="top center",
                textfont=dict(size=14, color='black', family='Arial Black'),
                name=f"C{label}: {nombre}",
                hovertemplate=f'<b>Corriente {label}</b><br>{nombre}<br>H: {h:.2f} kJ/mol<br>P: {p:.3f} MPa<extra></extra>',
                showlegend=True
            ))

        fig.update_layout(
            title="Diagrama Presión-Entalpía del Ciclo de Brayton",
            xaxis_title="Entalpía (kJ/mol)",
            yaxis_title="Presión (MPa)",
            yaxis_type="log",
            hovermode='closest',
            height=600,
            showlegend=False
        )

        st.plotly_chart(fig, use_container_width=True)

        # Información del ciclo
        with st.expander("ℹ️ Descripción del Ciclo"):
            st.markdown("""
            **Secuencia del Ciclo de Brayton:**

            **Ciclo de potencia (mezcla CO₂+H₂O):**
            - **C3** → Productos de combustión a alta P y T (entrada turbina)
            - **C4** → Gases expandidos a la salida de la turbina
            - **C5** → Gases enfriados después del recuperador (mezcla CO₂+H₂O)

            **Discontinuidad (línea roja punteada):**
            - **C5 → C7** → Separador de agua (cambio de composición, NO proceso termodinámico continuo)

            **Ciclo de recirculación (CO₂ puro):**
            - **C7** → CO₂ puro separado
            - **C8** → CO₂ comprimido (supercrítico)
            - **C9** → CO₂ enfriado
            - **C10** → CO₂ recirculado
            - **C11** → CO₂ precalentado que retorna a la cámara de combustión

            **Nota:** La línea roja punteada entre C5 y C7 representa el separador, que NO es un proceso
            termodinámico del ciclo Brayton, sino una separación física de componentes (extracción de H₂O líquido).
            """)
    else:
        st.info("Presiona el botón 'SIMULAR' para ver el diagrama P-H.")

# ============================================================================
# TAB 4: Diagrama T-S
# ============================================================================
if tab4.is_active:
    st.header("Diagrama Temperatura-Entropía")

    if simulacion_exitosa and simulador:
        # Corrientes del ciclo: 3 -> 4 -> 5 -> [separador] -> 7 -> 8 -> 9 -> 10 -> 11
        # NOTA: Eliminamos C6 porque el separador NO es un proceso termodinámico continuo
        corrientes_ciclo = [3, 4, 5, 7, 8, 9, 10, 11]

        S_values = []
        T_values = []
        labels = []
        nombres = []

        # Colores para el ciclo principal (mismo esquema que P-H)
        color_map = {
            3: '#d62728',   # Rojo - Productos combustión (entrada turbina)
            4: '#ff7f0e',   # Naranja - Salida turbina
            5: '#bcbd22',   # Amarillo-verde - Salida recuperador (mezcla CO2+H2O)
            7: '#9467bd',   # Púrpura - CO2 puro (después del separador)
            8: '#8c564b',   # Café - CO2 comprimido
            9: '#e377c2',   # Rosa - CO2 enfriado
            10: '#7f7f7f',  # Gris - CO2 recirculado
            11: '#2ca02c'   # Verde - CO2 precalentado
        }

        for i in corrientes_ciclo:
            if i in simulador.corrientes:
                c = simulador.corrientes[i]
                if c.s and c.T:
                    S_values.append(c.s / 1000)  # kJ/(mol·K)
                    T_values.append(c.T - 273.15)  # °C
                    labels.append(f"{i}")
                    nombres.append(c.nombre)

        # Crear gráfico con Plotly
        fig = go.Figure()

        # Línea conectando C3->C4->C5 (ciclo de potencia con mezcla)
        if len(S_values) >= 3:
            fig.add_trace(go.Scatter(
                x=S_values[0:3],  # C3, C4, C5
                y=T_values[0:3],
                mode='lines',
                line=dict(color='rgba(100, 100, 100, 0.5)', width=2, dash='solid'),
                showlegend=False,
                hoverinfo='skip'
            ))

        # Discontinuidad C5->C7 (SEPARADOR - cambio de composición)
        if len(S_values) >= 4:
            idx_c5 = 2  # C5 es el tercer elemento (índice 2)
            idx_c7 = 3  # C7 es el cuarto elemento (índice 3)
            fig.add_trace(go.Scatter(
                x=[S_values[idx_c5], S_values[idx_c7]],
                y=[T_values[idx_c5], T_values[idx_c7]],
                mode='lines',
                line=dict(color='rgba(255, 0, 0, 0.4)', width=2, dash='dot'),  # Rojo punteado
                showlegend=False,
                hoverinfo='skip'
            ))

        # Línea conectando C7->C8->C9->C10->C11 (ciclo de CO2 puro)
        if len(S_values) >= 4:
            fig.add_trace(go.Scatter(
                x=S_values[3:],  # C7, C8, C9, C10, C11
                y=T_values[3:],
                mode='lines',
                line=dict(color='rgba(100, 100, 100, 0.5)', width=2, dash='solid'),
                showlegend=False,
                hoverinfo='skip'
            ))

        # Línea cerrando el ciclo (11->3) - Verde punteado
        if len(S_values) >= 2:
            fig.add_trace(go.Scatter(
                x=[S_values[-1], S_values[0]],  # Del último (C11) al primero (C3)
                y=[T_values[-1], T_values[0]],
                mode='lines',
                line=dict(color='rgba(0, 128, 0, 0.6)', width=3, dash='dash'),  # Verde punteado
                showlegend=False,
                hoverinfo='skip'
            ))

        # Puntos individuales con colores y etiquetas
        for i, (s, t, label, nombre) in enumerate(zip(S_values, T_values, labels, nombres)):
            corriente_id = corrientes_ciclo[i]
            fig.add_trace(go.Scatter(
                x=[s],
                y=[t],
                mode='markers+text',
                marker=dict(size=16, color=color_map[corriente_id], line=dict(width=2, color='black')),
                text=f"<b>{label}</b>",
                textposition="top center",
                textfont=dict(size=14, color='black', family='Arial Black'),
                name=f"C{label}: {nombre}",
                hovertemplate=f'<b>Corriente {label}</b><br>{nombre}<br>S: {s:.4f} kJ/(mol·K)<br>T: {t:.1f} °C<extra></extra>',
                showlegend=True
            ))

        fig.update_layout(
            title="Diagrama Temperatura-Entropía del Ciclo de Brayton",
            xaxis_title="Entropía (kJ/mol·K)",
            yaxis_title="Temperatura (°C)",
            hovermode='closest',
            height=600,
            showlegend=False
        )

        st.plotly_chart(fig, use_container_width=True)

        # Información del ciclo
        with st.expander("ℹ️ Descripción del Ciclo"):
            st.markdown("""
            **Secuencia del Ciclo de Brayton:**

            **Ciclo de potencia (mezcla CO₂+H₂O):**
            - **C3** → Productos de combustión a alta P y T (entrada turbina)
            - **C4** → Gases expandidos a la salida de la turbina
            - **C5** → Gases enfriados después del recuperador (mezcla CO₂+H₂O)

            **Discontinuidad (línea roja punteada):**
            - **C5 → C7** → Separador de agua (cambio de composición, NO proceso termodinámico continuo)

            **Ciclo de recirculación (CO₂ puro):**
            - **C7** → CO₂ puro separado
            - **C8** → CO₂ comprimido (supercrítico)
            - **C9** → CO₂ enfriado
            - **C10** → CO₂ recirculado
            - **C11** → CO₂ precalentado que retorna a la cámara de combustión

            **Nota:** La línea roja punteada entre C5 y C7 representa el separador, que NO es un proceso
            termodinámico del ciclo Brayton, sino una separación física de componentes (extracción de H₂O líquido).
            """)
    else:
        st.info("Presiona el botón 'SIMULAR' para ver el diagrama T-S.")

# ============================================================================
# TAB 5: Análisis de Sensibilidad
# ============================================================================
if tab5.is_active:
    st.header("🔬 Análisis de Sensibilidad: Fracción de Recirculación")

    st.markdown("""
    Este análisis evalúa cómo la **fracción de recirculación de CO₂** afecta:
    - La **eficiencia global del ciclo** (incluye turbina, compresor CO₂, ASU y compresor combustible)
    - La **potencia generada** por la turbina
    - La **composición de los productos** (fracción de agua en C4)
    """)

    # Controles del análisis
    st.subheader("Configuración del Análisis")

    col1, col2, col3 = st.columns(3)

    with col1:
        combustible_analisis = st.selectbox(
            "Combustible a analizar",
            ["Gas Natural", "Gas de Síntesis", "Propano", "Etanol"],
            help="Combustible a usar en el análisis de sensibilidad"
        )

    with col2:
        f_min = st.number_input(
            "Fracción mínima (%)",
            min_value=0.0,
            max_value=99.0,
            value=0.0,
            step=5.0,
            help="Fracción de recirculación mínima (debe ser < fracción máxima)"
        )

        f_max = st.number_input(
            "Fracción máxima (%)",
            min_value=0.0,
            max_value=99.0,
            value=95.0,
            step=5.0,
            help="Fracción de recirculación máxima (recomendado: ≤97% para estabilidad)"
        )

    with col3:
        num_pasos = st.number_input(
            "Número de pasos",
            min_value=3,
            max_value=50,
            value=10,
            step=1,
            help="Cantidad de puntos a simular entre min y max"
        )

    # Validación de entradas
    if f_min >= f_max:
        st.error("⚠️ Error: La fracción mínima debe ser menor que la máxima")
        analisis_valido = False
    elif f_max >= 100:
        st.error("⚠️ Error: La fracción máxima debe ser menor a 100%")
        analisis_valido = False
    else:
        analisis_valido = True

    # Botón para ejecutar análisis
    if st.button("▶️ CORRER ANÁLISIS", type="primary", disabled=not analisis_valido):
        with st.spinner("Ejecutando análisis de sensibilidad... Esto puede tomar unos momentos."):

            # Importar numpy para el análisis
            import numpy as np
            import pandas as pd

            # Vector de fracciones a analizar
            fracciones = np.linspace(f_min/100, f_max/100, num_pasos)

            # Arrays para almacenar resultados
            resultados = {
                'f_recirculacion': [],
                'eta_CCS': [],
                'W_turb': [],
                'x_H2O_C4': [],
                'T_combustion': [],  # Temperatura de combustión calculada
                'corrientes': []  # Lista de diccionarios con todas las corrientes
            }

            # Barra de progreso
            progress_bar = st.progress(0)
            status_text = st.empty()

            # Ejecutar simulaciones
            for i, f_recirc in enumerate(fracciones):
                status_text.text(f"Simulando paso {i+1}/{num_pasos} (f = {f_recirc*100:.1f}%)")

                # Crear parámetros para esta simulación (copia de parametros original)
                # NOTA: P_recirculacion = P_combustion (restricción termodinámica)
                # NOTA: T_separador se calcula automáticamente basado en P_salida_turbina
                T_sep_auto = calcular_T_separador_automatica(P_salida_turbina * 1e6)

                params_analisis = {
                    'tipo_combustible': combustible_analisis,
                    # 'metodo_mezcla' eliminado: automático (HEOS para puros, PR para mezclas)
                    'T_ambiente': T_ambiente + 273.15,
                    'P_combustion': P_combustion * 1e6,  # Convertir de MPa a Pa
                    'T_combustion': T_combustion + 273.15,
                    'P_salida_turbina': P_salida_turbina * 1e6,  # Convertir de MPa a Pa
                    'eta_turbina': eta_turbina,
                    'eta_compresor_fuel': eta_compresor_fuel,
                    'eta_compresor_CO2': eta_compresor_CO2,
                    'efectividad_recuperador': efectividad_recuperador,
                    'flujo_combustible': flujo_combustible,
                    'fraccion_recirculacion': f_recirc,  # Variable del análisis
                    'P_recirculacion': P_combustion * 1e6,  # Convertir de MPa a Pa (igual a P_combustion)
                    'T_separador': T_sep_auto  # Calculada automáticamente
                }

                try:
                    # Crear y ejecutar simulador
                    sim = SimuladorBrayton(params_analisis)
                    exito = sim.simular()

                    if exito:
                        # Almacenar resultados
                        resultados['f_recirculacion'].append(f_recirc * 100)
                        resultados['eta_CCS'].append(sim.eta_CCS)
                        resultados['W_turb'].append(sim.W_turbina)  # Trabajo de turbina (no W_neto)

                        # Obtener temperatura de combustión
                        T_comb = sim.T_combustion_calculada - 273.15 if hasattr(sim, 'T_combustion_calculada') else None
                        resultados['T_combustion'].append(T_comb)

                        # Obtener fracción de agua en C4
                        if 4 in sim.corrientes and sim.corrientes[4].composicion:
                            x_H2O = sim.corrientes[4].composicion.get("H2O", 0) * 100
                        else:
                            x_H2O = 0
                        resultados['x_H2O_C4'].append(x_H2O)

                        # Guardar todas las corrientes para análisis detallado
                        corrientes_dict = {}
                        for num_corriente, corriente in sim.corrientes.items():
                            corrientes_dict[num_corriente] = {
                                'T': corriente.T - 273.15 if corriente.T else None,  # °C
                                'P': corriente.P / 1e5 if corriente.P else None,  # bar
                                'h': corriente.h / 1000 if corriente.h else None,  # kJ/mol
                                's': corriente.s if corriente.s else None,  # J/(mol·K)
                                'rho': corriente.rho if corriente.rho else None,  # kg/m³
                                'flujo': corriente.flujo_molar if corriente.flujo_molar else None,  # mol/s
                                'composicion': corriente.composicion.copy() if corriente.composicion else {}
                            }
                        resultados['corrientes'].append(corrientes_dict)
                    else:
                        # Simulación falló, usar None
                        resultados['f_recirculacion'].append(f_recirc * 100)
                        resultados['eta_CCS'].append(None)
                        resultados['W_turb'].append(None)
                        resultados['x_H2O_C4'].append(None)
                        resultados['T_combustion'].append(None)
                        resultados['corrientes'].append({})

                except Exception as e:
                    st.warning(f"⚠️ Error en simulación con f={f_recirc*100:.1f}%: {str(e)}")
                    resultados['f_recirculacion'].append(f_recirc * 100)
                    resultados['eta_CCS'].append(None)
                    resultados['W_turb'].append(None)
                    resultados['x_H2O_C4'].append(None)
                    resultados['T_combustion'].append(None)
                    resultados['corrientes'].append({})

                # Actualizar barra de progreso
                progress_bar.progress((i + 1) / num_pasos)

            status_text.text("✅ Análisis completado")
            progress_bar.empty()

            # Almacenar resultados en session_state
            st.session_state['resultados_sensibilidad'] = resultados
            st.session_state['combustible_analisis'] = combustible_analisis

    # Mostrar resultados si existen
    if 'resultados_sensibilidad' in st.session_state:
        resultados = st.session_state['resultados_sensibilidad']
        combustible_usado = st.session_state.get('combustible_analisis', 'N/A')

        st.success(f"✅ Resultados del análisis con **{combustible_usado}**")

        # ====================================================================
        # GRÁFICO: eta_CCS y W_turb vs x_H2O o f_recirculacion
        # ====================================================================
        st.subheader("📈 Gráfico de Sensibilidad")

        # Selector de tipo de gráfico
        tipo_grafico = st.radio(
            "Seleccionar eje horizontal:",
            options=["Fracción de H₂O en C4", "Fracción de Recirculación"],
            horizontal=True,
            help="Elige qué variable deseas ver en el eje X del gráfico"
        )

        # Crear figura con ejes secundarios
        from plotly.subplots import make_subplots

        # Determinar datos del eje X según selección
        if tipo_grafico == "Fracción de H₂O en C4":
            x_data = resultados['x_H2O_C4']
            x_title = "Fracción molar de H₂O en C4 (%)"
            x_range = [min(resultados['x_H2O_C4'])-2, max(resultados['x_H2O_C4'])+2]
            subtitle = "Eficiencia Global y Potencia de Turbina vs Fracción de Agua en C4"
        else:  # "Fracción de Recirculación"
            x_data = resultados['f_recirculacion']
            x_title = "Fracción de Recirculación (%)"
            x_range = [min(resultados['f_recirculacion'])-2, max(resultados['f_recirculacion'])+2]
            subtitle = "Eficiencia Global y Potencia de Turbina vs Fracción de Recirculación"

        fig = make_subplots(
            specs=[[{"secondary_y": True}]],
            subplot_titles=[subtitle]
        )

        # Agregar traza de eficiencia global (eje izquierdo)
        fig.add_trace(
            go.Scatter(
                x=x_data,
                y=resultados['eta_CCS'],
                name="η_Global (%)",
                mode='lines+markers',
                marker=dict(size=8, color='blue'),
                line=dict(width=2, color='blue')
            ),
            secondary_y=False
        )

        # Agregar traza de potencia turbina (eje derecho)
        fig.add_trace(
            go.Scatter(
                x=x_data,
                y=resultados['W_turb'],
                name="W_turbina (MW)",
                mode='lines+markers',
                marker=dict(size=8, color='red'),
                line=dict(width=2, color='red')
            ),
            secondary_y=True
        )

        # Configurar ejes
        fig.update_xaxes(
            title_text=x_title,
            range=x_range
        )
        fig.update_yaxes(title_text="Eficiencia Global (%)", secondary_y=False, color='blue')
        fig.update_yaxes(title_text="Potencia Turbina (MW)", secondary_y=True, color='red')

        fig.update_layout(
            height=500,
            hovermode='x unified',
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="right",
                x=1
            )
        )

        st.plotly_chart(fig, use_container_width=True)

        # Información adicional
        with st.expander("📊 Ver datos numéricos completos"):
            df_resumen = pd.DataFrame({
                'f_recirculación (%)': resultados['f_recirculacion'],
                'x_H₂O en C4 (%)': [f"{x:.2f}" if x is not None else "N/A" for x in resultados['x_H2O_C4']],
                'T_combustión (°C)': [f"{x:.1f}" if x is not None else "N/A" for x in resultados['T_combustion']],
                'η_Global (%)': [f"{x:.2f}" if x is not None else "N/A" for x in resultados['eta_CCS']],
                'W_turbina (MW)': [f"{x:.3f}" if x is not None else "N/A" for x in resultados['W_turb']]
            })
            st.dataframe(df_resumen, width="stretch")

        # ====================================================================
        # TABLAS POR CORRIENTE
        # ====================================================================
        st.subheader("📋 Análisis por Corriente")

        # Selector de corriente
        corrientes_disponibles = []
        if resultados['corrientes']:
            # Obtener corrientes del primer resultado exitoso
            for corrientes_dict in resultados['corrientes']:
                if corrientes_dict:
                    corrientes_disponibles = sorted(corrientes_dict.keys())
                    break

        if corrientes_disponibles:
            corriente_seleccionada = st.selectbox(
                "Seleccionar corriente",
                corrientes_disponibles,
                format_func=lambda x: f"Corriente {x}",
                help="Elige la corriente que deseas analizar en detalle"
            )

            # Construir tabla para la corriente seleccionada
            datos_corriente = {
                'f_recirculación (%)': resultados['f_recirculacion'],
                'T (°C)': [],
                'P (bar)': [],
                'h (kJ/mol)': [],
                's (J/mol·K)': [],
                'ρ (kg/m³)': [],
                'Flujo (mol/s)': []
            }

            # Agregar composición si existe
            composiciones = {}

            for corrientes_dict in resultados['corrientes']:
                if corriente_seleccionada in corrientes_dict:
                    c = corrientes_dict[corriente_seleccionada]
                    datos_corriente['T (°C)'].append(f"{c['T']:.2f}" if c['T'] is not None else "N/A")
                    datos_corriente['P (bar)'].append(f"{c['P']:.2f}" if c['P'] is not None else "N/A")
                    datos_corriente['h (kJ/mol)'].append(f"{c['h']:.2f}" if c['h'] is not None else "N/A")
                    datos_corriente['s (J/mol·K)'].append(f"{c['s']:.2f}" if c['s'] is not None else "N/A")
                    datos_corriente['ρ (kg/m³)'].append(f"{c['rho']:.2f}" if c['rho'] is not None else "N/A")
                    datos_corriente['Flujo (mol/s)'].append(f"{c['flujo']:.4f}" if c['flujo'] is not None else "N/A")

                    # Agregar composición
                    for componente, fraccion in c['composicion'].items():
                        if componente not in composiciones:
                            composiciones[componente] = []
                        composiciones[componente].append(f"{fraccion*100:.2f}%")
                else:
                    # Corriente no disponible en esta simulación
                    datos_corriente['T (°C)'].append("N/A")
                    datos_corriente['P (bar)'].append("N/A")
                    datos_corriente['h (kJ/mol)'].append("N/A")
                    datos_corriente['s (J/mol·K)'].append("N/A")
                    datos_corriente['ρ (kg/m³)'].append("N/A")
                    datos_corriente['Flujo (mol/s)'].append("N/A")

                    for componente in composiciones.keys():
                        composiciones[componente].append("N/A")

            # Agregar composiciones a los datos
            for componente, valores in composiciones.items():
                datos_corriente[f"x_{componente}"] = valores

            # Crear DataFrame y mostrar
            df_corriente = pd.DataFrame(datos_corriente)

            st.markdown(f"### Corriente {corriente_seleccionada}")
            st.dataframe(df_corriente, width="stretch")

            # Botón para descargar datos
            csv = df_corriente.to_csv(index=False).encode('utf-8')
            st.download_button(
                label=f"📥 Descargar datos de Corriente {corriente_seleccionada}",
                data=csv,
                file_name=f"corriente_{corriente_seleccionada}_sensibilidad.csv",
                mime="text/csv"
            )
        else:
            st.warning("⚠️ No hay datos de corrientes disponibles. Ejecuta el análisis primero.")

    else:
        st.info("👆 Configura los parámetros arriba y presiona **CORRER ANÁLISIS** para comenzar.")

# ============================================================================
# TAB 6: Optimización
# ============================================================================
if tab6.is_active:
    st.header("🎯 Optimización de Parámetros del Ciclo")
    st.markdown("""
    Esta herramienta optimiza automáticamente los parámetros del ciclo para **maximizar la eficiencia global** (incluye turbina, compresor CO₂, ASU y compresor combustible).

    **Restricción Termodinámica:**
    - ✅ **P_combustion = P_recirculacion ≥ 7.377 MPa** (Pc del CO₂)
    - El CO₂ se mantiene supercrítico en C8-C11 (compresión, enfriamiento y recirculación)
    - Después de la turbina (C4-C7), el CO₂ opera en condiciones subcríticas
    """)

    # Controles de optimización
    col_opt1, col_opt2 = st.columns(2)

    with col_opt1:
        combustible_nombre = st.selectbox(
            "Combustible para optimización",
            ["Gas Natural", "Gas de Síntesis", "Propano", "Etanol"],
            key="combustible_opt"
        )

        # Mapeo de nombres a tipos de combustible para el simulador
        mapeo_combustible = {
            "Gas Natural": "Gas Natural",
            "Gas de Síntesis": "Gas de Síntesis",
            "Propano": "Propano",
            "Etanol": "Etanol"
        }
        combustible_opt = mapeo_combustible[combustible_nombre]

    with col_opt2:
        max_evaluaciones_input = st.number_input(
            "Número máximo de evaluaciones",
            min_value=100,
            max_value=10000,
            value=1500,
            step=100,
            help="Número máximo de evaluaciones de la función objetivo que delimita el alcance de la optimización. " +
                 "Mayor valor = mejor exploración pero más tiempo. " +
                 "Rango sugerido: 500-1500 (rápido, ~2-5 min), 1500-3000 (balanceado, ~5-15 min), 3000-10000 (exhaustivo, ~15-45 min)",
            key="max_eval_opt"
        )

    # Mostrar información sobre los rangos de optimización
    with st.expander("ℹ️ Ver rangos de parámetros optimizables"):
        st.markdown("""
        **Parámetros que serán optimizados:**

        | Parámetro | Rango Mínimo | Rango Máximo | Notas |
        |-----------|--------------|--------------|-------|
        | **Presión de combustión / recirculación** | **8.1 MPa** (81 bar) | 30 MPa (300 bar) | Mín = 1.1×Pc(CO₂) para garantizar régimen supercrítico |
        | Presión salida turbina | 0.1 MPa (1 bar) | 0.5 MPa (5 bar) | Rango típico para ciclos Brayton |
        | Fracción de recirculación | 50% | **97%** | Límite realista (>97% inestable numéricamente) |
        | Flujo de combustible | 10 mol/s | 25 mol/s | Rango estándar para oxy-combustión |

        **Parámetros NO optimizables (fijos o automáticos):**
        - **T_ambiente**: Condición de frontera del sitio (definida en sidebar)
        - **T_separador**: Calculada automáticamente = T_sat_H₂O(P_salida_turbina) - 10°C
        - **T_combustión**: Calculada mediante balance energético (restricción < 1800°C aplicada)
        - **W_neto_mínimo**: Restricción de potencia mínima = 5.0 MW (evita soluciones no prácticas)

        ⚠️ **Nota**: P_combustion = P_recirculacion (restricción termodinámica). Mínimo = 1.1×Pc(CO₂) para régimen supercrítico estable

        **IMPORTANTE**: P_recirculacion y P_combustion son **iguales** por restricción termodinámica.
        El CO₂ recirculado no puede cambiar de presión en el recuperador (solo intercambia calor).

        **Física del proceso:**
        - **C3→C4 (Turbina)**: P cae de ~10-30 MPa a 0.1-1 MPa → CO₂ subcrítico ✅ Normal
        - **C4→C7 (Separador)**: Baja presión, CO₂ gas/subcrítico ✅ Esperado
        - **C7→C8 (Compresor)**: Comprime hasta P_combustion ≥ 7.377 MPa → CO₂ supercrítico ✅ Requerido
        - **C8-C11**: CO₂ permanece a P_combustion (supercrítico) para enfriamiento y recirculación
        """)

    # Botón para ejecutar optimización
    if st.button("🚀 EJECUTAR OPTIMIZACIÓN", key="btn_optimizar", type="primary"):
        # Contenedor para el mensaje de progreso (se puede limpiar después)
        info_container = st.empty()
        info_container.info("⏳ Optimización en progreso... Esto puede tomar varios minutos.")

        # Crear barra de progreso
        progress_bar = st.progress(0)
        status_text = st.empty()

        # Parámetros constantes (eficiencias de equipos)
        params_fijos = {
            'eta_turbina': eta_turbina,
            'eta_compresor_fuel': eta_compresor_fuel,
            'eta_compresor_CO2': eta_compresor_CO2,
            'efectividad_recuperador': efectividad_recuperador
            # 'metodo_mezcla' eliminado: automático (HEOS para puros, PR para mezclas)
        }

        # Constantes críticas de CO2
        Tc_CO2 = 304.13  # K (31°C)
        Pc_CO2 = 7.377e6  # Pa (73.77 bar)

        # NOTA IMPORTANTE: En un ciclo Brayton ABIERTO, el CO2 NO es supercrítico en toda la operación:
        # - En la TURBINA (C3→C4): P cae de ~100+ bar a ~1-10 bar → CO2 subcrítico es NORMAL
        # - En la salida y separador (C4-C7): P baja, CO2 subcrítico/gas
        # - En COMPRESIÓN para recirculación (C7→C8): Aquí SÍ queremos CO2 supercrítico
        #
        # Por lo tanto, las restricciones solo aplican a:
        # 1. P_recirculacion: debe ser > Pc para almacenar CO2 supercrítico
        # 2. Opcionalmente P_combustion si queremos trabajar en régimen supercrítico en cámara

        # Definir límites de parámetros optimizables
        # NOTA: P_recirculacion = P_combustion (restricción termodinámica)
        # NOTA: T_separador se calcula automáticamente (no optimizable)
        # NOTA: T_ambiente NO es optimizable (condición de frontera ambiental)
        # NOTA: T_combustion se calcula mediante balance energético (no optimizable)
        # [P_combustion, P_salida_turbina, f_recirculacion, flujo_combustible]
        bounds = [
            (Pc_CO2 * 1.1, 300e5),            # P_combustion: 8.1-300 bar (>Pc_CO2 para garantizar régimen supercrítico)
            (1.0e5, 5.0e5),                   # P_salida_turbina: 1-5 bar (rango típico)
            (0.50, 0.97),                     # f_recirculacion: 50-97% (límite realista, evita inestabilidad >97%)
            (10.0, 25.0),                     # flujo_combustible: 10-25 mol/s (rango estándar para Brayton con oxy-combustión)
        ]

        # Configuración del optimizador basado en número de evaluaciones
        popsize = 15
        max_evaluaciones = max_evaluaciones_input

        # Calcular número de generaciones aproximado (para maxiter)
        # Se usa un número alto para que el límite real sea el de evaluaciones
        max_generaciones = max(1, int(max_evaluaciones / popsize)) + 10  # +10 de margen

        # Contador de evaluaciones
        iter_count = [0]
        debe_detener = [False]

        # Función objetivo
        # Objetivo: Maximizar eta_CCS (eficiencia global con CCS)
        # Restricciones:
        #   1. T_combustion < 1800°C (límite realista: turbinas avanzadas ~1650°C, Allam ~1100-1600°C)
        #   2. CO2 en corrientes de recirculación debe estar supercrítico (P > 7.377 MPa)
        #   3. W_neto >= 5 MW (evitar soluciones de alta eficiencia pero baja potencia)
        def objetivo(x):
            try:
                iter_count[0] += 1

                # Detener si se alcanzó el límite de evaluaciones
                if iter_count[0] > max_evaluaciones:
                    debe_detener[0] = True
                    return 1e10  # Retornar penalización alta para forzar detención

                # Actualizar cada 10 evaluaciones para mejor feedback
                if iter_count[0] % 10 == 0 or iter_count[0] == 1:
                    progress = min(iter_count[0] / max_evaluaciones, 1.0)
                    progress_bar.progress(progress)
                    status_text.text(f"Evaluación {iter_count[0]}/{max_evaluaciones}")

                # Extraer parámetros (ahora son 4: incluimos flujo_combustible)
                P_comb, P_salida_turb, f_recir, flujo_comb = x

                # Validación básica: P_combustion > P_salida_turbina
                if P_comb <= P_salida_turb:
                    return 1e10  # Penalización alta

                # Calcular T_separador automáticamente basado en P_salida_turbina
                T_sep_auto = calcular_T_separador_automatica(P_salida_turb)

                # Crear diccionario de parámetros para el simulador
                # NOTA: P_recirculacion = P_combustion (restricción termodinámica)
                # NOTA: T_separador se calcula automáticamente (garantiza condensación completa H2O)
                # NOTA: T_ambiente se usa del sidebar (condición de frontera, no optimizable)
                # NOTA: T_combustion se calcula mediante balance energético (no es input, es output)
                params_sim = {
                    'tipo_combustible': combustible_opt,
                    'T_ambiente': T_ambiente + 273.15,  # Del sidebar, en K
                    'P_combustion': P_comb,
                    'T_combustion': 1600.0 + 273.15,  # Estimación inicial (será recalculada)
                    'P_salida_turbina': P_salida_turb,
                    'P_recirculacion': P_comb,  # Igual a P_combustion
                    'T_separador': T_sep_auto,  # Calculada automáticamente
                    'flujo_combustible': flujo_comb,  # Optimizable
                    'fraccion_recirculacion': f_recir,
                    **params_fijos
                }

                # Crear simulador con parámetros
                sim = SimuladorBrayton(params_sim)

                # Ejecutar simulación
                exito = sim.simular()

                if not exito:
                    return 1e10  # Penalización si la simulación falla

                # RESTRICCIÓN CRÍTICA: T_combustion debe ser < 1800°C (límite realista materiales)
                # Basado en: turbinas avanzadas ~1650°C, Allam cycle ~1100-1600°C
                T_combustion_calculada = sim.corrientes[3].T if 3 in sim.corrientes else 3000.0
                T_max_combustion = 1800.0 + 273.15  # 1800°C en K (límite conservador realista)
                penalizacion_temperatura = 0
                if T_combustion_calculada > T_max_combustion:
                    # Penalización severa si T > 1800°C
                    exceso_T = T_combustion_calculada - T_max_combustion
                    penalizacion_temperatura = exceso_T * 1e5  # Penalización muy alta

                # Validar que CO2 esté supercrítico SOLO en corrientes de alta presión
                # DIAGRAMA CORREGIDO:
                # C7, C8, C12: Baja presión (NO validar, subcrítico es normal)
                # C9: CO2 comprimido a P_recirculacion → DEBE ser supercrítico
                # C10, C11: CO2 en sistema de recirculación → DEBEN ser supercríticos
                corrientes_recirculación = [9, 10, 11]  # Solo corrientes de alta presión
                penalizacion_supercritico = 0

                for nombre_corr in corrientes_recirculación:
                    if nombre_corr in sim.corrientes:
                        corr = sim.corrientes[nombre_corr]
                        # Verificar si tiene CO2
                        if 'CO2' in corr.composicion and corr.composicion['CO2'] > 0.5:
                            # Penalizar si no está supercrítico (solo validar presión, T puede variar)
                            if corr.P <= Pc_CO2:
                                # Penalización proporcional a la desviación de presión
                                dP = max(0, Pc_CO2 - corr.P)
                                penalizacion_supercritico += (dP/Pc_CO2) * 1e6

                # Calcular eficiencias
                eta_ciclo = sim.eta_cycle if sim.eta_cycle is not None else 0
                eta_O2 = sim.eta_O2 if sim.eta_O2 is not None else 0
                eta_CCS = sim.eta_CCS if sim.eta_CCS is not None else 0

                # RESTRICCIÓN: Potencia neta mínima (evitar soluciones de alta eficiencia pero baja potencia)
                W_neto = sim.W_neto if hasattr(sim, 'W_neto') and sim.W_neto is not None else 0
                W_minimo = 5.0  # MW mínimo requerido (ajustable según aplicación)
                penalizacion_potencia = 0
                if W_neto < W_minimo:
                    # Penalización proporcional al déficit de potencia
                    deficit_potencia = W_minimo - W_neto
                    penalizacion_potencia = deficit_potencia * 1e4  # Penalización alta

                # Función objetivo: maximizar eta_Global (eta_CCS) (minimizar -eta_CCS)
                # Agregar penalizaciones por restricciones
                objetivo_val = -eta_CCS + penalizacion_temperatura + penalizacion_supercritico + penalizacion_potencia

                return objetivo_val

            except Exception as e:
                # Si hay error, penalizar
                return 1e10

        # Ejecutar optimización con differential_evolution
        try:
            resultado = differential_evolution(
                objetivo,
                bounds,
                maxiter=max_generaciones,  # Número de generaciones
                popsize=popsize,           # Tamaño de población por generación
                strategy='best1bin',
                seed=42,
                disp=False,
                polish=False,  # Desactivar polish para control de evaluaciones
                workers=1,
                updating='deferred',
                atol=0,        # Tolerancia absoluta (0 = no parar por convergencia)
                tol=0.001      # Tolerancia relativa para convergencia
            )

            progress_bar.progress(1.0)
            status_text.text(f"Optimización completada! ({iter_count[0]} evaluaciones totales)")

            if resultado.success or resultado.fun < 1e9:
                # Extraer parámetros óptimos (ahora son 4: incluimos flujo_combustible)
                x_opt = resultado.x
                P_comb_opt, P_salida_opt, f_recir_opt, flujo_comb_opt = x_opt

                # Calcular T_separador automáticamente basado en P_salida_opt
                T_sep_opt = calcular_T_separador_automatica(P_salida_opt)

                # Ejecutar simulación final con parámetros óptimos
                # P_recirculacion = P_combustion (restricción termodinámica)
                # T_separador se calcula automáticamente (garantiza condensación completa H2O)
                # T_ambiente se usa del sidebar (condición de frontera, no optimizable)
                # T_combustion se calcula mediante balance energético (output, no input)
                params_optimo = {
                    'tipo_combustible': combustible_opt,
                    'T_ambiente': T_ambiente + 273.15,  # Del sidebar, en K
                    'P_combustion': P_comb_opt,
                    'T_combustion': 1600.0 + 273.15,  # Estimación inicial (será recalculada)
                    'P_salida_turbina': P_salida_opt,
                    'P_recirculacion': P_comb_opt,  # Igual a P_combustion
                    'T_separador': T_sep_opt,  # Calculada automáticamente
                    'flujo_combustible': flujo_comb_opt,  # Optimizable
                    'fraccion_recirculacion': f_recir_opt,
                    **params_fijos
                }
                sim_optimo = SimuladorBrayton(params_optimo)
                sim_optimo.simular()

                # Extraer T_combustion calculada
                T_comb_opt = sim_optimo.T_combustion_calculada if hasattr(sim_optimo, 'T_combustion_calculada') else (1600.0 + 273.15)

                # Guardar en session state
                st.session_state['sim_optimo'] = sim_optimo
                st.session_state['params_optimos'] = x_opt

                # Limpiar mensaje de progreso, barra y status
                info_container.empty()
                progress_bar.empty()
                status_text.empty()

                st.success("✅ Optimización completada exitosamente!")

            else:
                # Limpiar mensaje de progreso, barra y status
                info_container.empty()
                progress_bar.empty()
                status_text.empty()

                st.error("❌ La optimización no convergió. Intenta con más evaluaciones o ajusta los rangos de parámetros.")
                st.write(f"Razón: {resultado.message}")

        except Exception as e:
            # Limpiar mensaje de progreso, barra y status
            info_container.empty()
            progress_bar.empty()
            status_text.empty()

            st.error(f"❌ Error durante la optimización: {str(e)}")
            import traceback
            st.code(traceback.format_exc())

    # Mostrar resultados guardados en session_state (si existen)
    if 'sim_optimo' in st.session_state and st.session_state['sim_optimo'] is not None:
        st.markdown("---")
        st.header("📊 Resultados de Optimización")

        sim_optimo = st.session_state['sim_optimo']
        x_opt = st.session_state['params_optimos']
        P_comb_opt, P_salida_opt, f_recir_opt, flujo_comb_opt = x_opt

        # Calcular T_separador
        T_sep_opt = calcular_T_separador_automatica(P_salida_opt)
        T_comb_opt = sim_optimo.T_combustion_calculada if hasattr(sim_optimo, 'T_combustion_calculada') else (1600.0 + 273.15)

        # TABLA 1: Resumen de parámetros óptimos
        st.subheader("📊 Tabla 1: Parámetros Óptimos del Ciclo")

        # Mostrar parámetros fijos (no optimizables)
        st.markdown(f"""
        **Parámetros fijos/calculados (no optimizables):**
        - **T_ambiente:** {T_ambiente:.1f} °C (condición del sitio)
        - **T_separador:** {T_sep_opt - 273.15:.1f} °C (calculada = T_sat_H₂O - 10°C)
        - **T_combustión:** {T_comb_opt - 273.15:.1f} °C (calculada mediante balance energético, restricción: < 2000°C)
        """)

        params_data = {
            'Parámetro': [
                'Presión de combustión / recirculación',
                'Presión salida turbina',
                'Fracción de recirculación',
                'Flujo de combustible'
            ],
            'Valor Óptimo': [
                f"{P_comb_opt/1e6:.2f}",
                f"{P_salida_opt/1e6:.3f}",
                f"{f_recir_opt*100:.2f}",
                f"{flujo_comb_opt:.2f}"
            ],
            'Unidades': ['MPa', 'MPa', '%', 'mol/s']
        }

        df_params = pd.DataFrame(params_data)
        st.dataframe(df_params, width="stretch")

        # Botón de descarga para tabla de parámetros
        csv_params = df_params.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 Descargar parámetros óptimos",
            data=csv_params,
            file_name=f"parametros_optimos.csv",
            mime="text/csv",
            key="download_params_persistente"
        )

        # TABLA 2: Resumen de propiedades termodinámicas de corrientes
        st.subheader("📋 Tabla 2: Propiedades Termodinámicas de las Corrientes")

        corrientes_data = []
        for nombre, corriente in sim_optimo.corrientes.items():
            # Calcular composición como string
            comp_str = ", ".join([f"{comp}: {frac:.3f}" for comp, frac in corriente.composicion.items()])

            corrientes_data.append({
                'Corriente': nombre,
                'T [°C]': f"{corriente.T - 273.15:.2f}" if corriente.T else "N/A",
                'P [MPa]': f"{corriente.P/1e6:.3f}" if corriente.P else "N/A",
                'h [kJ/mol]': f"{corriente.h/1000:.2f}" if corriente.h else "N/A",
                's [J/(mol·K)]': f"{corriente.s:.2f}" if corriente.s else "N/A",
                'ρ [kg/m³]': f"{corriente.rho:.2f}" if corriente.rho else "N/A",
                'Flujo [mol/s]': f"{corriente.flujo_molar:.2f}" if corriente.flujo_molar else "N/A",
                'Composición': comp_str
            })

        df_corrientes = pd.DataFrame(corrientes_data)
        st.dataframe(df_corrientes, width="stretch")

        # Botón de descarga para tabla de corrientes
        csv_corrientes = df_corrientes.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 Descargar propiedades de corrientes",
            data=csv_corrientes,
            file_name=f"corrientes_optimizadas.csv",
            mime="text/csv",
            key="download_corrientes_persistente"
        )

        # TABLA 3: Eficiencias
        st.subheader("⚡ Tabla 3: Eficiencias del Ciclo")

        col_eff1, col_eff2, col_eff3 = st.columns(3)

        with col_eff1:
            st.metric(
                label="η Ciclo",
                value=f"{sim_optimo.eta_cycle:.2f} %"
            )

        with col_eff2:
            st.metric(
                label="η con ASU",
                value=f"{sim_optimo.eta_O2:.2f} %"
            )

        with col_eff3:
            st.metric(
                label="η Global",
                value=f"{sim_optimo.eta_CCS:.2f} %"
            )

        # Mostrar trabajos y potencias
        st.markdown("---")
        st.subheader("🔧 Balance Energético")

        # Calcular trabajos de compresores desde las corrientes
        if 1 in sim_optimo.corrientes:
            n_combustible = sim_optimo.corrientes[1].flujo_molar
            T_amb_K = sim_optimo.params['T_ambiente']
            P_amb = 101325  # Pa

            # Crear combustible y corriente de entrada temporal
            combustible_temp = Combustible(sim_optimo.params['tipo_combustible'])
            corriente_0_temp = Corriente(
                "Entrada compresor fuel",
                T=T_amb_K,
                P=P_amb,
                composicion=combustible_temp.composicion,
                flujo_molar=n_combustible
            )
            corriente_0_temp.calcular_propiedades()
            h_0 = corriente_0_temp.h  # J/mol
            h_1 = sim_optimo.corrientes[1].h  # J/mol
            W_comp_fuel = (n_combustible * (h_1 - h_0)) / 1e6  # MW
        else:
            W_comp_fuel = 0.0

        # Trabajo compresor CO2 recirculación
        # Intentar obtener del atributo calculado
        if hasattr(sim_optimo, 'W_CO2comp_recirculacion'):
            W_comp_CO2_recirc = sim_optimo.W_CO2comp_recirculacion
        # Si no existe, calcular desde corrientes 8 y 9
        elif 8 in sim_optimo.corrientes and 9 in sim_optimo.corrientes:
            n_CO2_recirc = sim_optimo.corrientes[8].flujo_molar
            h_8 = sim_optimo.corrientes[8].h
            h_9 = sim_optimo.corrientes[9].h
            if h_8 and h_9:
                W_comp_CO2_recirc = (n_CO2_recirc * (h_9 - h_8)) / 1e6  # MW
            else:
                W_comp_CO2_recirc = 0.0
        else:
            W_comp_CO2_recirc = 0.0

        # Potencia ASU
        W_ASU = sim_optimo.W_ASU if hasattr(sim_optimo, 'W_ASU') else 0.0

        # W_turbina y W_neto ya están en MW
        col_work1, col_work2, col_work3, col_work4, col_work5 = st.columns(5)

        with col_work1:
            st.metric("W Turbina", f"{sim_optimo.W_turbina:.2f} MW")

        with col_work2:
            st.metric("W Comp. Combustible", f"{W_comp_fuel:.2f} MW")

        with col_work3:
            st.metric("W Comp. CO2 Recirc", f"{W_comp_CO2_recirc:.2f} MW",
                     help="Solo se comprime fracción recirculada (diagrama corregido)")

        with col_work4:
            st.metric("W ASU", f"{W_ASU:.2f} MW")

        with col_work5:
            st.metric("W Neto", f"{sim_optimo.W_neto:.2f} MW")

    else:
        st.info("👆 Configura el combustible y el número de iteraciones, luego presiona **EJECUTAR OPTIMIZACIÓN**.")

# ============================================================================
# Footer (removido - más limpio sin footer)
# ============================================================================
