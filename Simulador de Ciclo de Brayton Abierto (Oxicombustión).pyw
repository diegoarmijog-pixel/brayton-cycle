import os
import sys
import subprocess
import tkinter as tk

def iniciar():
    # Crear una pequeña ventana emergente (Splash Screen)
    root = tk.Tk()
    root.title("Cargando Simulador")
    
    # Quitar los bordes de la ventana para que parezca un cuadro de carga limpio
    root.overrideredirect(True)
    
    # Configurar tamaño y centrar en la pantalla
    ancho = 350
    alto = 100
    x = (root.winfo_screenwidth() // 2) - (ancho // 2)
    y = (root.winfo_screenheight() // 2) - (alto // 2)
    root.geometry(f"{ancho}x{alto}+{x}+{y}")
    
    # Añadir un marco y texto
    frame = tk.Frame(root, highlightbackground="#0068C9", highlightthickness=2)
    frame.pack(fill=tk.BOTH, expand=True)
    
    label = tk.Label(frame, text="🚀 Iniciando Simulador de Ciclo Brayton...\n\nPor favor, espera unos segundos.", font=("Arial", 11))
    label.pack(expand=True)

    # Obtener el directorio y lanzar Streamlit en segundo plano
    directorio = os.path.dirname(os.path.abspath(__file__))
    subprocess.Popen(
        [sys.executable, "-m", "streamlit", "run", "simulador_brayton.py"],
        cwd=directorio,
        creationflags=0x08000000
    )
    
    # Cerrar esta ventanita automáticamente después de 6 segundos (6000 ms)
    root.after(6000, root.destroy)
    root.mainloop()

if __name__ == '__main__':
    iniciar()