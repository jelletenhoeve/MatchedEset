###
### All generics for the MicroRNA project
###



## accessor functions
setGeneric('name',      def = function(object) {standardGeneric('name')})
setGeneric("name<-",    def = function(object, value) standardGeneric("name<-")) 
setGeneric('esetNames', def = function(object) {standardGeneric('esetNames')})
setGeneric('eset',      def = function(object, id) {standardGeneric('eset')})
setGeneric('thresholds',def = function(object) {standardGeneric('thresholds')})

setGeneric('subsetByFeatures',  def = function(object, indexList) {standardGeneric('subsetByFeatures')})
setGeneric('mapEset',    def = function(object, eset.id, mapTo, mapping, ...) {standardGeneric('mapEset')})
setGeneric('subset',    def = function(object, i, j, ...) {standardGeneric('subset')})



## interaction functions
setGeneric('interactionNames',  def = function(object) {standardGeneric('interactionNames')})
setGeneric('interactions',      def = function(object, id) {standardGeneric('interactions')})
setGeneric('interactions<-',    def = function(object, id, value, ...) {standardGeneric('interactions<-')})
setGeneric('interactionCounts', def = function(object, genesets, ...) {standardGeneric('interactionCounts')})

## correlation functions
setGeneric('correlationMatrix', def = function(object) {standardGeneric('correlationMatrix')})
setGeneric('correlations',      def = function(object, id, ...) {standardGeneric('correlations')})




## plot functions
setGeneric('correlationBoxplot', def = function(object, id, ...) {standardGeneric('correlationBoxplot')})
setGeneric('plotCorrelations',   def = function(object, ...) {standardGeneric('plotCorrelations')})
setGeneric('plotEcdfs',          def = function(object, ...) {standardGeneric('plotEcdfs')})

setGeneric('abundancePlot',      def = function(object, id, ...) {standardGeneric('abundancePlot')})
setGeneric('plotProfile',        def = function(object, ...) {standardGeneric('plotProfile')})
setGeneric('kmPlot',             def = function(object, ...) {standardGeneric('kmPlot')})
setGeneric('kmPlotPheno',        def = function(object, ...) {standardGeneric('kmPlotPheno')})
setGeneric('wilcoxTest',         def = function(object, id, ...) {standardGeneric('wilcoxTest')})

setGeneric('thresholdVariationPlot',    def = function(object, id, ...) {standardGeneric('thresholdVariationPlot')})


## global test functions
setGeneric('globalTest',         def = function(object, ...) {standardGeneric('globalTest')})
setGeneric('globalTestPheno',    def = function(object, ...) {standardGeneric('globalTestPheno')})

