#ifndef ZSW_CONFIG_H
#define ZSW_CONFIG_H

#ifdef WIN32
#ifdef _EXPORTING
#define ZSW_API __declspec(dllexport)
#else
#define ZSW_API __declspec(dllimport)
#endif
#else
#define ZSW_API
#endif

#define USING_DOUBLE_PRECISION

#ifdef USING_DOUBLE_PRECISION
namespace zsw {
typedef double Scalar;
}
#else
namespace zsw {
typedef float Scalar;
}
#endif

namespace zsw {
	enum COPY_TYPE {
		SHALLOW_COPY,
		DEEP_COPY
	};
}

#endif /* ZSW__CONFIG_H */
