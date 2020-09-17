#include "Force.h"

void Force::registerObject(Type object) {
	_object = object;
}

double Force::getEffect(double timelapse) {
	return _object.effect(pos);
}