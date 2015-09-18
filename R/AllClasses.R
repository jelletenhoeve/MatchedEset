###
### S4 class definitions
###



## MatchedExpressionSet

setClass(
	Class = 'MatchedExpressionSet',

	representation = representation(
		expressionSets        = 'list',
		interactionDataFrames = 'list',
		correlationMatrix     = 'matrix',
		name                  = 'character',
		thresholds            = 'numeric'
	)
)
