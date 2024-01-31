## Lets build our first SummarizedExperiment object
library("SummarizedExperiment")
## ?SummarizedExperiment

## De los ejemplos en la ayuda oficial

## Creamos los datos para nuestro objeto de tipo SummarizedExperiment
## para 200 genes a lo largo de 6 muestras
nrows <- 200
ncols <- 6
## Números al azar de cuentas
set.seed(20210223)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
## Información de nuestros genes
rowRanges <- GRanges(
  rep(c("chr1", "chr2"), c(50, 150)),
  IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
  strand = sample(c("+", "-"), 200, TRUE),
  feature_id = sprintf("ID%03d", 1:200)
)
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))
## Información de nuestras muestras
colData <- DataFrame(
  Treatment = rep(c("ChIP", "Input"), 3),
  row.names = LETTERS[1:6]
)
## Juntamos ahora toda la información en un solo objeto de R
rse <- SummarizedExperiment(
  assays = SimpleList(counts = counts),
  rowRanges = rowRanges,
  colData = colData
)

## Exploremos el objeto resultante
rse

## Nombres de tablas de cuentas que tenemos (RPKM, CPM, counts, logcounts, etc)
assayNames(rse)


## El inicio de nuestra tabla de cuentas
head(assay(rse))


#### Ejercicio ###
## Comando 1
rse[1:2, ]
#Selecciona primeros dos renglones (primeros 2 genes).
#subconjunto de rowranges
## Comando 2
rse[, c("A", "D", "F")]
#hace como un "subset" de todo. Subconjunto de muestras.
#No es subconjunto de los genes


#Explorar objeto
library(iSEE)
iSEE::iSEE(rse)

### Ejercicio###
## Descarguemos unos datos de spatialLIBD
sce_layer <- spatialLIBD::fetch_data("sce_layer")
sce_layer
## Revisemos el tamaño de este objeto
lobstr::obj_size(sce_layer)
iSEE::iSEE(sce_layer)

#1. Descarga un PDF que reproduzca la imagen del lado derecho de
#la siguiente diapositiva. Incluye ese PDF en tu
#repositorio de notas del curso.
##

#2. Explora en con un heatmap la expresión de los genes MOBP, MBP y PCP4.
#Si hacemos un clustering (agrupamos los genes), ¿cúales genes se parecen más?
## En el clustering los genes más parecidos son MBP y MOBP

#3. ¿En qué capas se expresan más los genes MOBP y MBP?
##WM
