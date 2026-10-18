# Linked PRIVATE by Heavy Feather's own targets, so the flags never reach a host.
add_library(hf_warnings INTERFACE)

if(HF_WARNINGS)
  if(MSVC)
    target_compile_options(hf_warnings INTERFACE /W4 /permissive-)
    if(HF_WARNINGS_AS_ERRORS)
      target_compile_options(hf_warnings INTERFACE /WX)
    endif()
  else()
    target_compile_options(hf_warnings INTERFACE
      -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wnon-virtual-dtor -Wold-style-cast -Woverloaded-virtual)
    if(HF_WARNINGS_AS_ERRORS)
      target_compile_options(hf_warnings INTERFACE -Werror)
    endif()
  endif()
endif()
