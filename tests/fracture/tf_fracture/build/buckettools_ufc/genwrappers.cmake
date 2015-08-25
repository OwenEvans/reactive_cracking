execute_process(
    COMMAND /usr/local/terraferma/master/fenics-tferma-master/petsc-maint/reldebug/bin/systemwrappers_from_options -l /home/oevans/repos/bitbucket/git/reactive_cracking/tests/fracture/tf_fracture/build/buckettools_ufc/cpp_filenames.txt.temp -- /home/oevans/repos/bitbucket/git/reactive_cracking/tests/fracture/tf_fracture/dg_elasticity.tfml
    WORKING_DIRECTORY /home/oevans/repos/bitbucket/git/reactive_cracking/tests/fracture/tf_fracture/build/buckettools_ufc
    RESULT_VARIABLE RETCODE
    )

if (NOT RETCODE EQUAL 0)
  message(FATAL_ERROR "Command returned ${RETCODE}")
endif()

execute_process(
    COMMAND /usr/bin/cmake -E compare_files /home/oevans/repos/bitbucket/git/reactive_cracking/tests/fracture/tf_fracture/build/buckettools_ufc/cpp_filenames.txt.temp /home/oevans/repos/bitbucket/git/reactive_cracking/tests/fracture/tf_fracture/build/buckettools_ufc/cpp_filenames.txt
    RESULT_VARIABLE CPP_FILENAMES_CHANGED
    OUTPUT_QUIET
    ERROR_QUIET
    )

if (CPP_FILENAMES_CHANGED)
  execute_process(
      COMMAND /usr/bin/cmake -E copy_if_different /home/oevans/repos/bitbucket/git/reactive_cracking/tests/fracture/tf_fracture/build/buckettools_ufc/cpp_filenames.txt.temp /home/oevans/repos/bitbucket/git/reactive_cracking/tests/fracture/tf_fracture/build/buckettools_ufc/cpp_filenames.txt
      WORKING_DIRECTORY /home/oevans/repos/bitbucket/git/reactive_cracking/tests/fracture/tf_fracture/build/buckettools_ufc
      )

  execute_process(
      COMMAND /usr/bin/cmake --build /home/oevans/repos/bitbucket/git/reactive_cracking/tests/fracture/tf_fracture/build --target rebuild_cache
      WORKING_DIRECTORY /home/oevans/repos/bitbucket/git/reactive_cracking/tests/fracture/tf_fracture/build/buckettools_ufc
      )
endif (CPP_FILENAMES_CHANGED)
