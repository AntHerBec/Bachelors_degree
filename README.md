There are 3 main categories into which these codes fall:
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
1) PARTICLE ANALYSIS -> test.py, Plotter.py and Selector.py
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
They complement each other, importing defined classes and functions. Its function is to filter the signal (interesting particle collision events that produce the wanted particles, in this case top antitop
production) by reducing the background (processes we are not interested in), accounting for the necessary corrections to the data, simulating the remaining backgroung events, obtaining the different 
uncertainties and calculating the cross section of the top antitop process.
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
2) INFERENTIAL STUDY -> Estudio_Accidentes_VariablesTemporales.R and its source data (TABLA_ACCIDENTES_20.xlsx)
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
The only code in R in the repository. Its objective is to study and try to conclude if the provided "temporal variables" (hour of the day, day of the week and month of the year). First, a simple 
descriptive analysis in order to understand the defined random variables, and an inferential analysis applying various tests to decide whether to reject or not a number of hypotheses, such as if the number 
of victims is the same given every possible value of the variables.
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
3) PHYSICS SIMULATIONS -> All the others
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Each code is related to a different physical system or a different approach to it:
 - Desintegraciones.py: A Monte Carlo simulation of time evolution of the number of radioactive nuclei and a comparison with the theoretical model that describes the decay.
 - Fotomultiplicador.py: A Monte Carlo simulation of the efficiency and performance of a chain of photomultipliers, supposing an electronic production per electron following a poisson distribution.
 - Onda_EM_1D.py/Onda_EM_2D.py: Animated simulation of the temporal evolution of an electromagnetic wave after a perturbation by numerically approximating Maxwell's equations.
 - Función_de_Onda.py: Animated simualtion of the temporal evolution of both the real and imaginary parts of a wave function inside a certain potential well by numerical approximation.
 - Gas_Encerrado.py: Animated simulation of the behaviour of an ideal gas inside an isolated container, and calculation various physical magnitudes such as temperature and pressure.
 - Pozo_Potencial.py: Numerical approximation of the stable solutions to the wave equation under a finite square potential well and a simplified atomic potential, getting the allowed energies in the process.
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
