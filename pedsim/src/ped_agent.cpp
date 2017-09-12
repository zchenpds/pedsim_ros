//
// pedsim - A microscopic pedestrian simulation system.
// Copyright (c) 2003 - 20012 by Christian Gloor
//

#include "ped_agent.h"
#include "ped_obstacle.h"
#include "ped_scene.h"
#include "ped_waypoint.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <fstream>
#include <iostream>

using namespace std;

default_random_engine generator;

// HRI parameters
double A_ir = 4;
double B_ir = 1;
double r_r = 0.5;

// anticipatory force parameters
static const float _k = 1.5*3;
static const float _m = 2.0;
static const float _t0 = 3.0;
static const float MAX_FORCE = 20;
static const float EPSILON = 0.0001;

/// Default Constructor
Ped::Tagent::Tagent()
{
    static int staticid = 0;
    id = staticid++;
    p.x = 0;
    p.y = 0;
    p.z = 0;
    v.x = 0;
    v.y = 0;
    v.z = 0;
    type = ADULT;
    scene = nullptr;
    teleop = false;
    double mean, var;
        std::cout << "Using default pedvel distribution.";
        mean = 1.25;
        var = 0.0;
        forceFactorDesired = 1.0;
        forceFactorSocial = 1.0;

    // assign random maximal speed in m/s
    normal_distribution<double> distribution(mean, var);
    vmax = distribution(generator);
#if 0
    forceFactorDesired = 5.0;
    forceFactorSocial = 2.1;
#endif
    forceFactorObstacle = 0.5755; // 10.0; // Aw
    forceSigmaObstacle = 8.9152; // 0.2; // Bw

    agentRadius = 0.3;
    relaxationTime = 0.5009;
}

/// Destructor
Ped::Tagent::~Tagent() {}

/// Assigns a Tscene to the agent. Tagent uses this to iterate over all
/// obstacles and other agents in a scene.
/// The scene will invoke this function when Tscene::addAgent() is called.
/// \warning Bad things will happen if the agent is not assigned to a scene. But
/// usually, Tscene takes care of that.
/// \param   *s A valid Tscene initialized earlier.
void Ped::Tagent::assignScene(Ped::Tscene* sceneIn) { scene = sceneIn; }

void Ped::Tagent::removeAgentFromNeighbors(const Ped::Tagent* agentIn)
{
    // search agent in neighbors, and remove him
    set<const Ped::Tagent*>::iterator foundNeighbor = neighbors.find(agentIn);
    if (foundNeighbor != neighbors.end())
        neighbors.erase(foundNeighbor);
}

/// Sets the maximum velocity of an agent (vmax). Even if pushed by other
/// agents, it will not move faster than this.
/// \param   pvmax The maximum velocity. In scene units per timestep, multiplied
/// by the simulation's precision h.
void Ped::Tagent::setVmax(double pvmax) { vmax = pvmax; }

/// Sets the agent's position. This, and other getters returning coordinates,
/// will eventually changed to returning a
/// Tvector.
/// \param   px Position x
/// \param   py Position y
/// \param   pz Position z
void Ped::Tagent::setPosition(double px, double py, double pz)
{
    p.x = px;
    p.y = py;
    p.z = pz;
}

/// Sets the factor by which the desired force is multiplied. Values between 0
/// and about 10 do make sense.
/// \param   f The factor
void Ped::Tagent::setForceFactorDesired(double f) { forceFactorDesired = f; }

/// Sets the factor by which the social force is multiplied. Values between 0
/// and about 10 do make sense.
/// \param   f The factor
void Ped::Tagent::setForceFactorSocial(double f) { forceFactorSocial = f; }

/// Sets the factor by which the obstacle force is multiplied. Values between 0
/// and about 10 do make sense.
/// \param   f The factor
void Ped::Tagent::setForceFactorObstacle(double f) { forceFactorObstacle = f; }

/// Calculates the force between this agent and the next assigned waypoint.
/// If the waypoint has been reached, the next waypoint in the list will be
/// selected.
/// \return  Tvector: the calculated force
Ped::Tvector Ped::Tagent::desiredForce()
{
    // get destination
    Twaypoint* waypoint = getCurrentWaypoint();

    // if there is no destination, don't move
    if (waypoint == NULL) {
        desiredDirection = Ped::Tvector();
        Tvector antiMove = -v / relaxationTime;
        return antiMove;
    }

    // compute force
    Tvector force = waypoint->getForce(*this, &desiredDirection);

    return force;
}

/// Calculates the social force between this agent and all the other agents
/// belonging to the same scene.
/// It iterates over all agents inside the scene, has therefore the complexity
/// O(N^2). A better
/// agent storing structure in Tscene would fix this. But for small (less than
/// 10000 agents) scenarios, this is just
/// fine.
/// \return  Tvector: the calculated force
Ped::Tvector Ped::Tagent::socialForce() const
{
#define USE_CUSTOM_HRI 2

    // define relative importance of position vs velocity vector
    // (set according to Moussaid-Helbing 2009)
    const double lambdaImportance = 2.0;

    // define speed interaction
    // (set according to Moussaid-Helbing 2009)
    const double gamma = 0.35;

    // define speed interaction
    // (set according to Moussaid-Helbing 2009)
    const double n = 2;

    // define angular interaction
    // (set according to Moussaid-Helbing 2009)
    const double n_prime = 3;

#if 0
    // HRI parameters
    const double A = 23;
    const double B = 0.1;
    const double r = 0.6;
#endif
    //cout << " A = " << A;
    //cout << " B = " << B;
    //cout << " r = " << r << endl;

    Tvector force;
    for (const Ped::Tagent* other : neighbors) {
        // don't compute social force to yourself
        if (other->id == id)
            continue;

        // compute difference between both agents' positions
        Tvector diff = other->p - p;
        Tvector diffDirection = diff.normalized();

        // compute difference between both agents' velocity vectors
        // Note: the agent-other-order changed here
        Tvector velDiff = v - other->v;

#if USE_CUSTOM_HRI==1
        if (other->getType() == ROBOT)
        {
            Tvector forceHRI = -A_ir * exp(((r_r + agentRadius) - diff.length())/B_ir) * diffDirection;
            force += forceHRI;

        }
        else
        {
         double const_D0 = 0.31;
         double const_D1 = 0.45;
         double const_ij = const_D1/diff.length(); 
         double var_ij = pow(const_ij, 2);
         double lambda = 0.3387; // Chao: 0.25
         double A = 25; //0.4072; // Aij
         double B = 0.08; //0.1959; // Bij
         double Theta_ij = lambda + (1-lambda)*0.5*(1 + Tvector::dotProduct(v.normalized(), diff.normalized()));
         // 
         double gab_ij;
         Tvector tab_ij;
         tab_ij.x = diffDirection.y;
         tab_ij.y = - diffDirection.x;         
         if (diff.length() > 2*agentRadius)
         {
             gab_ij = 0;
         }
         else
         {
             gab_ij = 2*agentRadius - diff.length();   
         }   

         //Tvector forcePeds = -exp(-diff.length() / const_D0 + var_ij) * diffDirection * Theta_ij + 0 * gab_ij * Tvector::dotProduct(-velDiff, tab_ij) * tab_ij;
         Tvector forcePeds = -A * exp( (2*agentRadius - diff.length()) / B ) * diffDirection * Theta_ij;
         if (2*agentRadius > diff.length())
             forcePeds += -1.5e4 * (2*agentRadius > diff.length()) * diffDirection;
         force += forcePeds; 
        }
#elif USE_CUSTOM_HRI==2


	if (other->getType() == ROBOT)
        {
	    ////////////////////////////////////////////////////////////////////////////////////////
	
	    const float distanceSq = (other->p - p).lengthSquared();
	    float radiusSq = (r_r + agentRadius) * (r_r + agentRadius);
	    if (distanceSq != radiusSq)
	    {	
		// if agents are actually colliding use their separation distance 
		if (distanceSq < radiusSq)
		    radiusSq  = 0.99 * distanceSq;
					
		const Tvector w_ = other->p - p;
		const Tvector v_ = v - other->v;
		const float a = Tvector::dotProduct(v_, v_);
		const float b = Tvector::dotProduct(w_, v_);	
		const float c = Tvector::dotProduct(w_, w_) - radiusSq;
		float discr = b*b - a*c;
		if (discr > .0f && (a<-EPSILON || a>EPSILON))
		{
		    discr = sqrt(discr);
		    const float t = (b - discr) / a;
		    if (t>0){
			Tvector forceHRI = -_k*exp(-t/_t0)*(v_ - (b*v_ - a*w_)/discr)/(a*pow(t,_m))*(_m/t + 1/_t0);
			if (forceHRI.length() > MAX_FORCE)
			    forceHRI = forceHRI / forceHRI.length() * MAX_FORCE;
			force += forceHRI;
		    }
		}
	    }


	    ////////////////////////////////////////////////////////////////////////////////////////
            //Tvector forceHRI = -A_ir * exp(((r_r + agentRadius) - diff.length())/B_ir) * diffDirection;
            //force += forceHRI;

        }
        else
        {
         double const_D0 = 0.31;
         double const_D1 = 0.45;
         double const_ij = const_D1/diff.length(); 
         double var_ij = pow(const_ij, 2);
         double lambda = 0.3387; // Chao: 0.25
         double A = 25; //0.4072; // Aij
         double B = 0.08; //0.1959; // Bij
         double Theta_ij = lambda + (1-lambda)*0.5*(1 + Tvector::dotProduct(v.normalized(), diff.normalized()));
         // 
         double gab_ij;
         Tvector tab_ij;
         tab_ij.x = diffDirection.y;
         tab_ij.y = - diffDirection.x;         
         if (diff.length() > 2*agentRadius)
         {
             gab_ij = 0;
         }
         else
         {
             gab_ij = 2*agentRadius - diff.length();   
         }   

         //Tvector forcePeds = -exp(-diff.length() / const_D0 + var_ij) * diffDirection * Theta_ij + 0 * gab_ij * Tvector::dotProduct(-velDiff, tab_ij) * tab_ij;
         Tvector forcePeds = -A * exp( (2*agentRadius - diff.length()) / B ) * diffDirection * Theta_ij;
         if (2*agentRadius > diff.length())
             forcePeds += -1.5e4 * (2*agentRadius > diff.length()) * diffDirection;
         force += forcePeds; 
        }
#else
    // compute interaction direction t_ij
    Tvector interactionVector = lambdaImportance * velDiff + diffDirection;
    double interactionLength = interactionVector.length();
    Tvector interactionDirection = interactionVector / interactionLength;

    // compute angle theta (between interaction and position difference vector)
    Ped::Tangle theta = interactionDirection.angleTo(diffDirection);

    // compute model parameter B = gamma * ||D||
    double B = gamma * interactionLength;

    double thetaRad = theta.toRadian();
    double forceVelocityAmount = -exp(-diff.length() / B - (n_prime * B * thetaRad) * (n_prime * B * thetaRad));
    double forceAngleAmount = -theta.sign() * exp(-diff.length() / B - (n * B * thetaRad) * (n * B * thetaRad));

    Tvector forceVelocity = forceVelocityAmount * interactionDirection;
    Tvector forceAngle = forceAngleAmount * interactionDirection.leftNormalVector();

    force += forceVelocity + forceAngle;
#endif

    }
    return force;


}


#if USE_CUSTOM_HRI==2
/// Calculates the force between this agent and the nearest obstacle in this
/// scene.
/// Iterates over all obstacles == O(N).
/// \return  Tvector: the calculated force
Ped::Tvector Ped::Tagent::obstacleForce() const
{

    //////////////////////////////////////////////////////////
    // obstacle which is closest only
//    Ped::Tvector minDiff;
    float effectiveObstacleDistance = 2;
    const float _INFTY = 1e10;
    Ped::Tvector forceObstacle;
    for (const Tobstacle* obstacle : scene->obstacles) 
    {
        //const LineObstacle* obstacle = simEngine->getObstacle(i);
        Ped::Tvector closestPoint = obstacle->closestPoint(p);
        Ped::Tvector n_w = closestPoint - p;
	float d_w = n_w.lengthSquared();
	

	// Agent is moving away from obstacle, already colliding or obstacle too far away
	if (Ped::Tvector::dotProduct(v, n_w) < 0 || d_w == agentRadius * agentRadius || d_w > effectiveObstacleDistance *effectiveObstacleDistance)
	    continue;
        // correct the radius, if the agent is already colliding	
	const float radius = d_w < (agentRadius*agentRadius)? sqrt(d_w):agentRadius; 

	const float a= Ped::Tvector::dotProduct(v, v);
	bool discCollision = false, segmentCollision = false;
	float t_min = _INFTY;
	        
	float c, b, discr; 
	float b_temp, discr_temp, c_temp, D_temp;
	Ped::Tvector w_temp, w, o1_temp, o2_temp, o_temp, o, w_o;

	// time-to-collision with disc_1 of the capped rectangle (capsule)
	w_temp = obstacle->getStartPoint() - p;
	b_temp = Ped::Tvector::dotProduct(w_temp, v);	
	c_temp = Ped::Tvector::dotProduct(w_temp, w_temp) - (radius*radius);
	discr_temp = b_temp*b_temp - a*c_temp;
	if (discr_temp > .0f && (a<-EPSILON || a>EPSILON))
	{
	    discr_temp = sqrt(discr_temp);
	    const float t = (b_temp - discr_temp) / a;
	    if (t>0){
		t_min = t;
		b = b_temp;
		discr = discr_temp;
		w = w_temp;
		c = c_temp;
		discCollision = true;
	    }
	}
					
	// time-to-collision with disc_2 of the capsule
	w_temp = obstacle->getEndPoint() - p;
	b_temp = Ped::Tvector::dotProduct(w_temp, v);	
	c_temp = Ped::Tvector::dotProduct(w_temp, w_temp) - (radius*radius);
	discr_temp = b_temp*b_temp - a*c_temp;
	if (discr_temp > .0f && (a<-EPSILON || a>EPSILON))
	{
	    discr_temp = sqrt(discr_temp);
	    const float t = (b_temp - discr_temp) / a;
	    if (t>0 && t < t_min){
		t_min = t;
		b = b_temp;
		discr = discr_temp;
		w = w_temp;
		c = c_temp;
		discCollision = true;
	    }
	}

	// time-to-collision with segment_1 of the capsule
	Ped::Tvector p1 = obstacle->getStartPoint();
	Ped::Tvector p2 = obstacle->getEndPoint();
	
	o1_temp = p1 + radius * (p2 - p1).leftNormalVector();
	o2_temp = p2 + radius * (p2 - p1).leftNormalVector();
	o_temp = o2_temp - o1_temp;

	D_temp = Ped::Tvector::crossProduct(v, o_temp).z;
	if (D_temp != 0)   
	{
	    float inverseDet = 1.0f/D_temp;	
	    float t = Ped::Tvector::crossProduct(o_temp, p - o1_temp).z * inverseDet;
	    float s = Ped::Tvector::crossProduct(v, p - o1_temp).z * inverseDet;
	    if (t>0 && s>=0 && s <=1 && t<t_min)
	    {
		t_min = t;
		o = o_temp;
		w_o = p - o1_temp;
		discCollision = false;
		segmentCollision = true;
	    }
	}

	// time-to-collision with segment_2 of the capsule
	o1_temp = p1 + radius * (p2 - p1).rightNormalVector();
	o2_temp = p2 + radius * (p2 - p1).rightNormalVector();
	o_temp = o2_temp - o1_temp;

	D_temp = Ped::Tvector::crossProduct(v, o_temp).z;
	if (D_temp != 0)   
	{
	    float inverseDet = 1.0f/D_temp;	
	    float t = Ped::Tvector::crossProduct(o_temp, p - o1_temp).z * inverseDet;
	    float s = Ped::Tvector::crossProduct(v, p - o1_temp).z * inverseDet;
	    if (t>0 && s>=0 && s <=1 && t<t_min)
	    {
		t_min = t;
		o = o_temp;
		w_o = p - o1_temp;
		discCollision = false;
		segmentCollision = true;
	    }
	}

	if (discCollision)
	{
	    forceObstacle += -_k*exp(-t_min/_t0)*(v - (b*v - a*w)/discr)/(a*pow(t_min,_m))*(_m/t_min + 1/_t0);		
	}
	else if (segmentCollision)
	{
	    forceObstacle += _k*exp(-t_min/_t0)/(pow(t_min,_m)*Ped::Tvector::crossProduct(v,o)).z*(_m/t_min + 1/_t0)*Ped::Tvector(-o.y, o.x);
	}
        
    }
    if (forceObstacle.length() > MAX_FORCE)
	forceObstacle = forceObstacle / forceObstacle.length() * MAX_FORCE;

    // physical obstacle force
    // obstacle which is closest only
    Ped::Tvector minDiff;
    double minDistanceSquared = INFINITY;

    for (const Tobstacle* obstacle : scene->obstacles) {
        Ped::Tvector closestPoint = obstacle->closestPoint(p);
        Ped::Tvector diff = p - closestPoint;
        double distanceSquared = diff.lengthSquared(); // use squared distance to
        // avoid computing square
        // root
        if (distanceSquared < minDistanceSquared) {
            minDistanceSquared = distanceSquared;
            minDiff = diff;
        }
    }

    double distance = sqrt(minDistanceSquared) - agentRadius;
    double forceAmount = exp(-distance / forceSigmaObstacle);
    if (agentRadius > sqrt(minDistanceSquared))
             forceAmount += 1.5e2 * (agentRadius - sqrt(minDistanceSquared));
    forceObstacle += forceAmount * minDiff.normalized();

    if (forceObstacle.length() > MAX_FORCE)
	forceObstacle = forceObstacle / forceObstacle.length() * MAX_FORCE;
    return forceObstacle;
    //////////////////////////////////////////////////////////
}



#else
/// Calculates the force between this agent and the nearest obstacle in this
/// scene.
/// Iterates over all obstacles == O(N).
/// \return  Tvector: the calculated force
Ped::Tvector Ped::Tagent::obstacleForce() const
{
    // obstacle which is closest only
    Ped::Tvector minDiff;
    double minDistanceSquared = INFINITY;

    for (const Tobstacle* obstacle : scene->obstacles) {
        Ped::Tvector closestPoint = obstacle->closestPoint(p);
        Ped::Tvector diff = p - closestPoint;
        double distanceSquared = diff.lengthSquared(); // use squared distance to
        // avoid computing square
        // root
        if (distanceSquared < minDistanceSquared) {
            minDistanceSquared = distanceSquared;
            minDiff = diff;
        }
    }

    double distance = sqrt(minDistanceSquared) - agentRadius;
    double forceAmount = exp(-distance / forceSigmaObstacle);
    if (agentRadius > sqrt(minDistanceSquared))
             forceAmount += 1.5e2 * (agentRadius - sqrt(minDistanceSquared));
    return forceAmount * minDiff.normalized();
}
#endif


/// myForce() is a method that returns an "empty" force (all components set to
/// 0).
/// This method can be overridden in order to define own forces.
/// It is called in move() in addition to the other default forces.
/// \return  Tvector: the calculated force
/// \param   e is a vector defining the direction in which the agent wants to
/// walk to.
Ped::Tvector Ped::Tagent::myForce(Ped::Tvector e) const
{
    return Ped::Tvector();
}

void Ped::Tagent::computeForces()
{
    // update neighbors
    // NOTE - have a config value for the neighbor range
    const double neighborhoodRange = 10.0;
    neighbors = scene->getNeighbors(p.x, p.y, neighborhoodRange);

    // update forces
    desiredforce = desiredForce();
    if (forceFactorSocial > 0)
        socialforce = socialForce();
    if (forceFactorObstacle > 0)
        obstacleforce = obstacleForce();
    myforce = myForce(desiredDirection);
}

/// Does the agent dynamics stuff. Calls the methods to calculate the individual
/// forces, adds them
/// to get the total force affecting the agent. This will then be translated
/// into a velocity difference,
/// which is applied to the agents velocity, and then to its position.
/// \param   stepSizeIn This tells the simulation how far the agent should
/// proceed
void Ped::Tagent::move(double stepSizeIn)
{
    // sum of all forces --> acceleration
    a = forceFactorDesired * desiredforce + forceFactorSocial * socialforce + forceFactorObstacle * obstacleforce + myforce;
    /*           std::cout << "A_desired = " << forceFactorDesired << "\n";
               std::cout << "A =  " << forceFactorSocial << "\n";
               std::cout << "Aw = " << forceFactorObstacle << "\n Bw = " << forceSigmaObstacle * obstacleforce.length() << endl;*/
    // calculate the new velocity
    if (getTeleop() == false) {
        v = v + stepSizeIn * a;
    }

    // don't exceed maximal speed
    double speed = v.length();
    if (speed > vmax)
        v = v.normalized() * vmax;

    // internal position update = actual move
    p += stepSizeIn * v;

    // notice scene of movement
    scene->moveAgent(this);
}
