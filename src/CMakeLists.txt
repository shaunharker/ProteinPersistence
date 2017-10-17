# add the executable

set ( LIBS      ${LIBS}
                ${Boost_LIBRARIES}
                ${VTK_LIBRARIES}
                sqlite3 )

add_executable ( protein2diagram pers_weighted_alpha_3d.cpp )
target_link_libraries ( protein2diagram ${LIBS} )

install(TARGETS protein2diagram
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
