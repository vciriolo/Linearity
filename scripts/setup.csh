if (!(-d bin/)) then
mkdir lib
endif
if (!(-d lib/)) then
mkdir lib
endif
if (!(-d obj/)) then
mkdir lib
endif

if (${?LINEARITY}) then
echo "already set"
else
setenv THISDIR `pwd`

setenv LD_LIBRARY_PATH ${THISDIR}/../CommonUtils/lib:${LD_LIBRARY_PATH}
#setenv DYLD_LIBRARY_PATH ${THISDIR}/../CommonUtils/lib:${DYLD_LIBRARY_PATH}
setenv PATH ${THISDIR}/../CommonUtils/bin:${PATH}

setenv LD_LIBRARY_PATH ${THISDIR}/lib:${LD_LIBRARY_PATH}
#setenv DYLD_LIBRARY_PATH ${THISDIR}/lib:${DYLD_LIBRARY_PATH}
setenv PATH ${THISDIR}/bin:${PATH}

setenv COMMONUTILS ${THISDIR}/../CommonUtils
setenv COMMONUTILSINCLUDE ${THISDIR}/../CommonUtils/interface
setenv COMMONUTILSLIB ${THISDIR}/../CommonUtils/lib

setenv LINEARITY ${THISDIR}/
setenv LINEARITYINCLUDE ${THISDIR}/interface
setenv LINEARITYLIB ${THISDIR}/lib
endif
