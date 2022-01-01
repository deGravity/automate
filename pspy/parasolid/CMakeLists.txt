cmake_minimum_required(VERSION 3.8)

add_library(parasolid STATIC "frustrum.cpp" "parasolid.cpp")

IF (WIN32)
	IF (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/lib/pskernel_archive_win_x64.lib)
		execute_process(
			COMMAND cmake -E tar xf ${CMAKE_CURRENT_SOURCE_DIR}/lib/pskernel_archive_win_x64.zip
			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib
		)
	ENDIF()
	target_link_libraries(parasolid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/pskernel_archive_win_x64.lib)
ELSEIF(APPLE)
	set(THREADS_PREFER_PTHREAD_FLAG ON)
	find_package(Threads REQUIRED)
	IF (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/lib/pskernel_archive_intel_macos.lib)
		execute_process(
			COMMAND cmake -E tar xf ${CMAKE_CURRENT_SOURCE_DIR}/lib/pskernel_archive_intel_macos.zip
			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib
		)
	ENDIF()
	target_link_libraries(parasolid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/pskernel_archive_intel_macos.lib Threads::Threads)
ELSE()
	set(THREADS_PREFER_PTHREAD_FLAG ON)
	find_package(Threads REQUIRED)
	IF (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/lib/pskernel_archive_linux_x64.lib)
		execute_process(
			COMMAND cmake -E tar xf ${CMAKE_CURRENT_SOURCE_DIR}/lib/pskernel_archive_linux_x64.zip
			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib
		)
	ENDIF()
	target_link_libraries(parasolid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/pskernel_archive_linux_x64.lib Threads::Threads)
ENDIF()

target_include_directories(parasolid INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
