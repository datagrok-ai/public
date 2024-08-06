cmake_minimum_required(VERSION 3.10)

project(XGBoostAPI)

# Specify the C++ standard
add_definitions(-std=c++17)
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(xgboost)

if (DEFINED EMSCRIPTEN)
	add_executable(XGBoostAPI xgboost-api.cpp)

        target_link_libraries(XGBoostAPI xgboost)

	set_target_properties(XGBoostAPI PROPERTIES COMPILE_FLAGS "-O3")
	set_target_properties(XGBoostAPI PROPERTIES LINK_FLAGS    "-O3 -s TOTAL_MEMORY=268435456 -s WASM=1 -s ALLOW_MEMORY_GROWTH=1 -s MODULARIZE=1 -s EXPORT_NAME=\"XGBoost\" -s EXPORTED_FUNCTIONS=[\"_train\",\"_predict\",\"_malloc\",\"_free\"] -s EXPORTED_RUNTIME_METHODS=[\"cwrap\",\"ccall\"] -sENVIRONMENT=web,worker")
else()
	add_library(XGBoostApi SHARED xgboost-api.cpp)
endif()



