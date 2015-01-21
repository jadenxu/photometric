SET(GRAPH_CUT_LIBS
    block.h
    energy.h
    graph.h
    graph.cpp
    LinkedBlockList.h
    LinkedBlockList.cpp
    GCoptimization.h
    GCoptimization.cpp
    maxflow.cpp)
ADD_LIBRARY(graphCutCore ${GRAPH_CUT_LIBS})
