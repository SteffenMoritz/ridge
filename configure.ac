AC_INIT(ridge, version-2.9, steffen.moritz10@gmail.com, ridge_2.9.tar.gz, http://cran.r-project.org/web/packages/ridge)
#include <stdlib.h>

AC_CONFIG_HEADER([src/config.h])


dnl Set the variable HAVE_GSL to TRUE
HAVE_GSL=TRUE

# Checks for common programs using default macros
AC_PROG_CC

## Use gsl-config to find arguments for compiler and linker flags
##
## Check for non-standard programs: gsl-config(1)
AC_PATH_PROG([GSL_CONFIG], [gsl-config])
## If gsl-config was found, let's use it
if test "${GSL_CONFIG}" != ""; then
    # Use gsl-config for header and linker arguments
    GSL_CFLAGS=`${GSL_CONFIG} --cflags`
    GSL_LIBS=`${GSL_CONFIG} --libs`
else
    AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi




dnl Substitute HAVE_GSL in the functions R/l*Genotypes*.R.in 
AC_SUBST(HAVE_GSL)
AC_CONFIG_FILES([R/l*Genotypes*.R])
AC_OUTPUT


dnl Now substitute these variables in src/Makevars.in to create src/Makevars
AC_SUBST(GSL_CFLAGS) 
AC_SUBST(GSL_LIBS)
AC_OUTPUT(src/Makevars)

