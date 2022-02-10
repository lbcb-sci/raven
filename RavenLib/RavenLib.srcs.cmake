
set(SOURCES "")

set(SOURCES
        ${CMAKE_CURRENT_LIST_DIR}/Source/pile.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Source/pile.hpp
        ${CMAKE_CURRENT_LIST_DIR}/Source/Graph/Graph.hpp
        ${CMAKE_CURRENT_LIST_DIR}/Source/Graph/Graph.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Source/Graph/OverlapUtils.hpp
        ${CMAKE_CURRENT_LIST_DIR}/Source/Graph/GraphAssemble.hpp
        ${CMAKE_CURRENT_LIST_DIR}/Source/Graph/GraphAssemble.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Source/Graph/GraphPolish.hpp
        ${CMAKE_CURRENT_LIST_DIR}/Source/Graph/GraphPolish.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Source/Graph/GraphShared.hpp
        ${CMAKE_CURRENT_LIST_DIR}/Source/Graph/GraphShared.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Source/Graph/GraphConstruct.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Source/Graph/GraphConstruct.hpp
        ${CMAKE_CURRENT_LIST_DIR}/Source/Graph/Serialization/GraphBinarySerialization.hpp
        ${CMAKE_CURRENT_LIST_DIR}/Source/Graph/Serialization/GraphBinarySerialization.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Source/Graph/Serialization/GraphOutputs.hpp
        )

source_group("Sources" FILES ${SOURCES})
list(APPEND SOURCES ${Sources})
