// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdOdIobjdImydict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "../interface/ConfigFileLine.h"
#include "../interface/ConfigParser.h"
#include "../interface/GetCategories.h"
#include "../interface/GetExtraSmearing.h"
#include "../interface/GetScaleCorrection.h"
#include "../interface/LHEReader.h"
#include "../interface/PUReweighting.h"
#include "../interface/ParseIJazZFile.h"
#include "../interface/ScaleEstimators.h"
#include "../interface/ZPeakFitUtils.h"
#include "../interface/geometryUtils.h"
#include "../interface/histoFunc.h"
#include "../interface/ntpleUtils.h"
#include "../interface/setTDRStyle.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *vectorlEfloatgR_Dictionary();
   static void vectorlEfloatgR_TClassManip(TClass*);
   static void *new_vectorlEfloatgR(void *p = 0);
   static void *newArray_vectorlEfloatgR(Long_t size, void *p);
   static void delete_vectorlEfloatgR(void *p);
   static void deleteArray_vectorlEfloatgR(void *p);
   static void destruct_vectorlEfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<float>*)
   {
      vector<float> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<float>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<float>", -2, "vector", 214,
                  typeid(vector<float>), DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<float>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<float>*)0x0)->GetClass();
      vectorlEfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEfloatgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<float> : new vector<float>;
   }
   static void *newArray_vectorlEfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<float>[nElements] : new vector<float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEfloatgR(void *p) {
      delete ((vector<float>*)p);
   }
   static void deleteArray_vectorlEfloatgR(void *p) {
      delete [] ((vector<float>*)p);
   }
   static void destruct_vectorlEfloatgR(void *p) {
      typedef vector<float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<float>

namespace {
  void TriggerDictionaryInitialization_mydict_Impl() {
    static const char* headers[] = {
"../interface/ConfigFileLine.h",
"../interface/ConfigParser.h",
"../interface/GetCategories.h",
"../interface/GetExtraSmearing.h",
"../interface/GetScaleCorrection.h",
"../interface/LHEReader.h",
"../interface/PUReweighting.h",
"../interface/ParseIJazZFile.h",
"../interface/ScaleEstimators.h",
"../interface/ZPeakFitUtils.h",
"../interface/geometryUtils.h",
"../interface/histoFunc.h",
"../interface/ntpleUtils.h",
"../interface/setTDRStyle.h",
0
    };
    static const char* includePaths[] = {
"../interface/",
"/afs/cern.ch/user/v/vciriolo/work/private/ECALlinearity/CMSSW_7_6_3/src/linearity_studies/Linearity/../CommonUtils/interface",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc493/cms/cmssw/CMSSW_7_6_3/external/slc6_amd64_gcc493/bin/../../../../../../lcg/root/6.02.12-kpegke4/include",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc493/lcg/root/6.02.12-kpegke4/include",
"/afs/cern.ch/work/v/vciriolo/private/ECALlinearity/CMSSW_7_6_3/src/linearity_studies/Linearity/src/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$string")))  allocator;
}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "../interface/ConfigFileLine.h"
#include "../interface/ConfigParser.h"
#include "../interface/GetCategories.h"
#include "../interface/GetExtraSmearing.h"
#include "../interface/GetScaleCorrection.h"
#include "../interface/LHEReader.h"
#include "../interface/PUReweighting.h"
#include "../interface/ParseIJazZFile.h"
#include "../interface/ScaleEstimators.h"
#include "../interface/ZPeakFitUtils.h"
#include "../interface/geometryUtils.h"
#include "../interface/histoFunc.h"
#include "../interface/ntpleUtils.h"
#include "../interface/setTDRStyle.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("mydict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_mydict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_mydict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_mydict() {
  TriggerDictionaryInitialization_mydict_Impl();
}
