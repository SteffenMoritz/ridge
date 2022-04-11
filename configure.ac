
AC_INIT([ridge],[version-3.3],[steffen.moritz10@gmail.com],[ridge_3.3.tar.gz],[http://cran.r-project.org/web/packages/ridge])

AC_CONFIG_HEADERS([src/config.h])


: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

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
    HAVE_GSL=TRUE
else
    AC_MSG_WARN([gsl-config not found, is GSL installed?])
    AC_MSG_WARN([ridge will be installed but some functions will be unavailable])
    HAVE_GSL=FALSE
fi

if test "${HAVE_GSL}" = TRUE; then
dnl Define HAVE_GSL_H in src/config.h
    AC_CHECK_HEADERS([gsl/gsl_version.h],
    AC_DEFINE(HAVE_GSL_HEADER,[], [Description]), [HAVE_GSL=FALSE])
    if test "${HAVE_GSL}" = FALSE; then
       AC_MSG_WARN([gsl headers not found, perhaps check the path?])
       AC_MSG_WARN([ridge will be installed, but some functions will be unavailable])
    fi
fi


dnl Substitute HAVE_GSL in the functions R/l*Genotypes*.R.in 
AC_SUBST(HAVE_GSL)
AC_CONFIG_FILES([R/l*Genotypes*.R])




dnl Now substitute these variables in src/Makevars.in to create src/Makevars
if test "${HAVE_GSL}" = TRUE; then
      AC_SUBST(GSL_CFLAGS) 
      AC_SUBST(GSL_LIBS)
      AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
fi