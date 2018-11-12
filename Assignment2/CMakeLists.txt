cmake_minimum_required(VERSION 3.8)
project(Assignment2)

set(CMAKE_CXX_STANDARD 11)

find_package( OpenCV REQUIRED )

include_directories( ${OpenCV_INCLUDE_DIRS} )

set(SOURCE main.cpp)

add_executable(${PROJECT_NAME} ${SOURCE})

target_link_libraries( ${PROJECT_NAME} ${OpenCV_LIBS} )