AC_INIT(ridge, version-2.8, steffen.moritz10@gmail.com, ridge_2.8.tar.gz, http://cran.r-project.org/web/packages/ridge)
AC_CONFIG_HEADER([src/config.h])
#include <stdlib.h>

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`


dnl Set the variable HAVE_GSL to TRUE
HAVE_GSL=TRUE
dnl Following lines kindly supplied by Dirk Eddelbuettel
dnl Use gsl-config to find arguments for compiler and linker flags
dnl Check for non-standard programs: gsl-config(1)
AC_PATH_PROG([GSL_CONFIG], [gsl-config])
dnl If gsl-config was found, let's use it
if test "${GSL_CONFIG}" != ""; then
dnl Use gsl-config for header and linker arguments of GSL 
    GSL_CFLAGS=`${GSL_CONFIG} --cflags`
    GSL_LIBS=`${GSL_CONFIG} --libs`
    CFLAGS="$CFLAGS $GSL_CFLAGS"
    CPPFLAGS="$CPPFLAGS $GSL_CPPFLAGS"
dnl else set HAVE_GSL to FALSE and print warnings
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

if test "${HAVE_GSL}" = TRUE; then
   AC_MSG_CHECKING([if GSL version >= 1.14])
   AC_RUN_IFELSE([AC_LANG_SOURCE([[
   #include <string.h>
   #include <gsl/gsl_version.h>
   int main() {
       #ifdef GSL_VERSION
       int major, minor;
       char *gslv = GSL_VERSION;
       if ((sscanf(gslv, "%d.%d", &major, &minor)) != 2) {
       	  exit (1);
       }
       exit (minor < 14);
       #else
       exit(1);
       #endif
       }
       ]])],
       [gsl_version_ok=yes],
       [gsl_version_ok=no],
       [gsl_version_ok=yes])
   if test "${gsl_version_ok}" = no; then
      	AC_MSG_RESULT([no])
	AC_MSG_WARN([gsl < 1.14 was found and will not be used])
	AC_MSG_WARN([ridge will be installed but some functions will be unavailable])
	AC_MSG_WARN([to enable GSL functions please upgrade GSL to >=1.14 and reinstall ridge])
	HAVE_GSL=FALSE	    
   else
	AC_MSG_RESULT([yes])
   fi dnl ends if test ${gsl_version_ok}
fi dnl ends if $HAVE_GSL

dnl Substitute HAVE_GSL in the functions R/l*Genotypes*.R.in 
AC_SUBST(HAVE_GSL)
AC_CONFIG_FILES([R/l*Genotypes*.R])
AC_OUTPUT

if test "${HAVE_GSL}" = TRUE; then
dnl Check whether libopenblas is available
AC_CHECK_LIB(openblas, $sgemm, [HAVE_OPENBLAS=TRUE;HAVE_OPENBLAS=FALSE])
dnl If openblas is available
if test "${HAVE_OPENBLAS}" = TRUE; then
   GSL_CBLAS_LIB=-lopenblas
   GSL_CFLAGS=`${GSL_CONFIG} --cflags`
   GSL_LIBS="`${GSL_CONFIG} --libs-without-cblas` -lopenblas"
fi
fi

dnl Now substitute these variables in src/Makevars.in to create src/Makevars
if test "${HAVE_GSL}" = TRUE; then
      AC_SUBST(GSL_CFLAGS) 
      AC_SUBST(GSL_LIBS)
      AC_OUTPUT(src/Makevars)
fi
	
