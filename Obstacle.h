#pragma once

#include "IObstacle.h"
#include "InteractingBody.h"

namespace MPCD {
	class Obstacle :
		public IObstacle,
		public InteractingBody
	{


	};

}

