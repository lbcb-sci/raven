
set(SOURCES "")

set(SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/src/assemble.cc 
  ${CMAKE_CURRENT_LIST_DIR}/src/binary.cc 
  ${CMAKE_CURRENT_LIST_DIR}/src/common.cc 
  ${CMAKE_CURRENT_LIST_DIR}/src/construct.cc 
  ${CMAKE_CURRENT_LIST_DIR}/src/graph_repr.cc
  ${CMAKE_CURRENT_LIST_DIR}/src/graph.cc 
  ${CMAKE_CURRENT_LIST_DIR}/src/io.cc
  ${CMAKE_CURRENT_LIST_DIR}/src/overlap_utils.cc
  ${CMAKE_CURRENT_LIST_DIR}/src/pile.cc 
  ${CMAKE_CURRENT_LIST_DIR}/src/polish.cc 
)

source_group("Sources" FILES ${SOURCES})
list(APPEND SOURCES ${Sources})
