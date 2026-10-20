# A stand-in host: built on its own against an INSTALLED Heavy Feather, never as part of the
# main tree. It proves the package works the way a game would use it:
#
#   cmake --install build/release --prefix /tmp/hf
#   cmake -S tests/package -B build/package -DCMAKE_PREFIX_PATH=/tmp/hf
#   cmake --build build/package && ./build/package/hf-consumer

cmake_minimum_required(VERSION 3.25)
project(HeavyFeatherConsumer LANGUAGES CXX)

find_package(HeavyFeather 0.3 REQUIRED)

add_executable(hf-consumer main.cpp)
target_link_libraries(hf-consumer PRIVATE HeavyFeather::HeavyFeather)
