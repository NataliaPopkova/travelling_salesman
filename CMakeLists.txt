cmake_minimum_required(VERSION 3.16)	 

project(travelling_salesman)			# Название проекта

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_EXE main.cpp)		# Установка переменной со списком исходников для исполняемого файла

set(SOURCE_LIB map.cpp)			# Тоже самое, но для библиотеки

add_library(map STATIC ${SOURCE_LIB})	# Создание статической библиотеки с именем map

add_executable(main ${SOURCE_EXE})	# Создает исполняемый файл с именем main

target_link_libraries(main map)		# Линковка программы с библиотекой