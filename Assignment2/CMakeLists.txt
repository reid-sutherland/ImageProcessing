cmake_minimum_required(VERSION 3.12)
project(Assignment2)

set(CMAKE_CXX_STANDARD 11)

find_package( OpenCV REQUIRED )

include_directories( ${OpenCV_INCLUDE_DIRS} )

add_executable(Assignment2 main.cpp)

target_link_libraries( ${PROJECT_NAME} ${OpenCV_LIBS} )