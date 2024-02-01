## ?model.matrix
#datos, y, x, x+1
mat <- with(trees, model.matrix(log(Volume) ~ log(Height) + log(Girth)))
mat

colnames(mat)

class(mat)
dim(mat)
dim(trees)
head(trees)
head(mat)

summary(lm(log(Volume) ~ log(Height) + log(Girth), data = trees))

#Bi es el cambio en el volumen de los árboles ajustados por
#el cambio en el log de la altura y de girth.

#Se puede usar para datos humanos. Ej. expresión de un
#gen ajustado a la edad.

#Guardar variable que nos interesa:
#Estimado/SD = t

#### Ejemplo ExploreModelMatrix ###
## Datos de ejemplo
(sampleData <- data.frame(
  genotype = rep(c("A", "B"), each = 4),
  treatment = rep(c("ctrl", "trt"), 4)
))

## Creemos las imágenes usando ExploreModelMatrix
vd <- ExploreModelMatrix::VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ genotype + treatment,
  textSizeFitted = 4
)

## Veamos las imágenes
cowplot::plot_grid(plotlist = vd$plotlist)

## Usaremos shiny otra ves
library(ExploreModelMatrix)

app <- ExploreModelMatrix(
  sampleData = sampleData,
  designFormula = ~ genotype + treatment
)
if (interactive()) shiny::runApp(app)

#Interfaz nos permite interpretar los coeficientes.
#Columna ii (con tratamineto), i (control)
#En este caso: Tratamiento trt: se cancelan al restar las celdas
#Columna ii es cuando el tratamiento trt está presente
#Pero con genotipo constante.
#E[Y|trat=trt,gen]-E[Y|treat=ctrl,gen]
#Muestras con tratamiento trt contra tramiento ctrl
#Con genotipo constante
#Valor Esperado = media = promedio

#Función ModelMatrix:
#Indicator dummy: 1 (no referencia), 0 (variable de referencia)


### Ejemplo 2 ####
#Datos con id de individuos que fueron muestreados
#antes y después del tratamineto
(sampleData <- data.frame(
  Response = rep(c("Resistant", "Sensitive"), c(12, 18)),
  Patient = factor(rep(c(1:6, 8, 11:18), each = 2)),
  Treatment = factor(rep(c("pre","post"), 15)),
  ind.n = factor(rep(c(1:6, 2, 5:12), each = 2))))

app <- ExploreModelMatrix(
  sampleData = sampleData,
  designFormula = ~ Response + Response:ind.n + Response:Treatment
)
if (interactive()) shiny::runApp(app)

#columna de referencia: post tratamiento. ii(pre), post(i)
#ii-i
#E[Y|treatment=pre,response=sentitive]-
#E[Y|treatment=post,response=sensitive]
#cambio promendio en la expresión, con individuos constantes,
#para individuos sensibles.

### Ejemplo 3 ###
(sampleData = data.frame(
  condition = factor(rep(c("ctrl_minus", "ctrl_plus",
                           "ko_minus", "ko_plus"), 3)),
  batch = factor(rep(1:6, each = 2))))
app <- ExploreModelMatrix(sampleData = sampleData,
                          designFormula = ~ batch + condition)
if (interactive()) shiny::runApp(app)

# ii-i con batch constante.
#E[Y|condition=ctrl_+, batch=1,3,5]-E[Y|condition=ctrl_-, batch=1,3.5]
#cambio en respuesta al tratamiento con batch constante.

### Ejercicio ###

#Interpreta ResponseResistant.Treatmentpre del ejercicio 2.
#Puede ser útil tomar un screenshot
#(captura de pantalla) y anotar la con líneas de colores. Si haces eso,
#puedes incluir la imagen en tus notas.
##E[y|treatment=pre, resp=resistant, indin] -
##E[y|treatment=post,resp=resistant,ind]

#¿Por qué es clave el 0 al inicio de la fórmula en el ejercicio 3?
## Es necesaria para descartar la columna intercept.
##Para no tener el promedio de todos los ratones. Y solo tener
##al batch1.
