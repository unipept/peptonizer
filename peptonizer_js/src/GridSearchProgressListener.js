// Interface that should be implemented by whoever wants to be notified about progress updates.
class GridSearchProgressListener {
    gridUpdated(currentAlpha, currentBeta, currentPrior) {}
    graphsUpdated(currentGraph, totalGraphs) {}
    maxResidualUpdated(maxResidual, tolerance) {}
    iterationsUpdated(currentIteration, totalIterations) {}
}
