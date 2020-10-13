
include_directories(includes/math)
include_directories(includes/general)

add_local_sources(sof impnse_fFilt.c)
add_local_sources(sof impnse_stdmlib.c)
add_local_sources(sof impnse_main.c)
add_local_sources(sof impnse.c)
add_local_sources(sof impnse_generic.c)


set(SOURCE_FILES 
   impnse_fFilt.c
   impnse_main.c
   impnse_stdmlib.c
   impnse.c
   impnse_generic.c)

add_executable(SOURCE_FILES ${SOURCE_FILES})
target_link_libraries(SOURCE_FILES PRIVATE -lm)
