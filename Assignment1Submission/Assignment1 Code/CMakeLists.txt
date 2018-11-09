cmake_minimum_required(VERSION 2.8)
project( assignment1 )

find_package( OpenCV REQUIRED )

set(CMAKE_CXX_STANDARD 11)

include_directories( ${OpenCV_INCLUDE_DIRS} )

set( SOURCES
        main.cpp
        )

add_executable( ${PROJECT_NAME} ${SOURCES} )
target_link_libraries( ${PROJECT_NAME} ${OpenCV_LIBS})