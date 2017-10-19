# add the executable

set ( LIBS      ${LIBS}
                ${VTK_LIBRARIES} )

add_executable ( pers_weighted_alpha_3d pers_weighted_alpha_3d.cpp )
target_link_libraries ( pers_weighted_alpha_3d ${LIBS} )

install(TARGETS pers_weighted_alpha_3d
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
