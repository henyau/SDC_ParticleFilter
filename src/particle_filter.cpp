/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 *  Implementation: Henry Yau
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;



void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Sets the number of particles. Initialize all particles to first position (based on estimates of 
	// x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Adds random Gaussian noise to each particle.
	//Consult particle_filter.h for more information about this method (and others in this file).
	
	default_random_engine generator;
	normal_distribution<double> distribution_x(x, std[0]);
	normal_distribution<double> distribution_y(y, std[1]);
	normal_distribution<double> distribution_theta(theta, std[2]);

	num_particles = 100; // test for now.

	Particle tempParticle;

	for (int i = 0; i < num_particles; i++)
	{
		tempParticle.x = distribution_x(generator);
		tempParticle.y = distribution_y(generator);
		tempParticle.theta = distribution_theta(generator);
		tempParticle.weight = 1.0f;
		tempParticle.id = i; // is this used?

		particles.push_back(tempParticle);
		weights.push_back(1.0f);
		
	}
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Adds measurements to each particle and add random Gaussian noise.
	
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine generator;
	normal_distribution<double> distribution_x(0, std_pos[0]);
	normal_distribution<double> distribution_y(0, std_pos[1]);
	normal_distribution<double> distribution_theta(0, std_pos[2]);

	double theta_0;
	for (int i = 0; i < num_particles; i++)
	{
		theta_0 = particles[i].theta;
		//update x,y
		if (fabs(yaw_rate) > 1E-8)
		{
			particles[i].x += (velocity / yaw_rate)*(sin(theta_0 + yaw_rate*delta_t) - sin(theta_0));
			particles[i].y += (velocity / yaw_rate)*(-cos(theta_0 + yaw_rate*delta_t) + cos(theta_0));
			
		}
		else
		{
			particles[i].x += velocity*cos(theta_0)*delta_t;
			particles[i].y += velocity*sin(theta_0)*delta_t;
		}

		//update yaw
		particles[i].theta += yaw_rate*delta_t;

		//add Gaussian noise
		particles[i].x += distribution_x(generator);
		particles[i].y += distribution_y(generator);
		particles[i].theta += distribution_theta(generator);
	}


}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	//Finds the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark	
//called for each particle, when done observations should contain actual landmark locations

	vector<LandmarkObs>::iterator predIter;
	vector<LandmarkObs>::iterator obsIter;
	double distObser;
	double distObserMin;
	LandmarkObs tempObs;
	
	vector<LandmarkObs> LandmarkClosest;

	//find nearest neighbor in prediction list for each observation
	for(obsIter = observations.begin(); obsIter != observations.end(); ++obsIter)
	{
		distObserMin = 1E20;
		for (predIter = predicted.begin(); predIter != predicted.end(); ++predIter)
		{
			distObser = dist(obsIter->x, obsIter->y, predIter->x, predIter->y);
			
			if (distObserMin > distObser)
			{
				distObserMin = distObser;
				// this is the closest observed landmark to the particle
				//obsIter->id = predIter->id;
				tempObs = *predIter;
			}
		}
		// tempObs is the closest observed landmark to the particle
		LandmarkClosest.push_back(tempObs);
		// we still need to retain a copy of the original transformed observation to compute the weights
	}
	observations = LandmarkClosest;
}
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// Updates the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	double dist_Particle;
	LandmarkObs tempObs;


	// predict measurements to landmarks for each particle
	for (int i = 0; i < num_particles; i++)
	{
//		vector<LandmarkObs>::iterator obsIter;
//		vector<LandmarkObs>::iterator predIter;

		std::vector<LandmarkObs> predVec;
		std::vector<LandmarkObs> obVec;
		std::vector<LandmarkObs> ob_w_LM_Vec; //nearest landmark to each observation

		//iterate through the map landmarks (in global)
		//for (obMapIter = map_landmarks.landmark_list.begin(); obMapIter != map_landmarks.landmark_list.end(); ++obMapIter)
		for(auto obMapIter: map_landmarks.landmark_list)
		{
			//distance from particle to observation			
			dist_Particle = dist(obMapIter.x_f, obMapIter.y_f, particles[i].x, particles[i].y);
			if (dist_Particle < sensor_range) // for every map location within the range, add to predVec
			{
				//if within sensor range, add to list of potential landmarks
				// this is in global coords
				tempObs.x = obMapIter.x_f;
				tempObs.y = obMapIter.y_f;
				tempObs.id = obMapIter.id_i;
				predVec.push_back(tempObs);
			}
		}
		
		//iterate through all the observations and transform to get a global coord
		//for (obsIter = observations.begin(); obsIter != observations.end(); ++obsIter)
		for (auto obsIter: observations)
		{
			//want the observation (from particle) in global coords
			homotransform(particles[i].x, particles[i].y, particles[i].theta, obsIter.x, obsIter.y, tempObs.x, tempObs.y);
			tempObs.id = obsIter.id;
			obVec.push_back(tempObs);
		}
		
		//get the closest landmark location from map (within range) for each observation
		ob_w_LM_Vec = obVec;
		dataAssociation(predVec, ob_w_LM_Vec);

		//update weight using the multivariate gaussian probability distribution

		double weight_accum = 1;
		double del_x, del_y;

		for (int i = 0; i<obVec.size(); i++)
		{
			del_x = obVec[i].x - ob_w_LM_Vec[i].x;
			del_y = obVec[i].y - ob_w_LM_Vec[i].y;
					
			weight_accum *= exp(-(del_x*del_x / (2 * std_landmark[0] * std_landmark[0]) + del_y*del_y/ (2 * std_landmark[1] * std_landmark[1]))) / (2 * M_PI*std_landmark[0] * std_landmark[1]);
				
		}
		//update the weight for each particle
		particles[i].weight = weight_accum;
		weights[i] = weight_accum; //keep seperate copy of weights to use in weighted rand dist
	}
	//update weights based on associations

}

void ParticleFilter::resample() {
	// Resamples particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	//
/*	
	default_random_engine generator;
	//create a new set of particles, select based on weights

	vector<Particle> newpList;
	int index = (int)floor(rand()*num_particles);
	double wmax = 0;
	for (int i = 0; i < num_particles; i++)
	{
		if (wmax < particles[i].weight)
			wmax = particles[i].weight;
	} 
	float beta = 0.0f;

	for (int i = 0; i < num_particles; i++)
	{
		beta += rand()*wmax * 2;
		do
		{
			beta -= particles[index].weight;
			index += 1;
			index = index%num_particles;
		} while (particles[index].weight < beta);
		newpList.push_back(particles[index]);
	}
	
	particles = newpList;
*/
	//std::random_device rd;
	//std::mt19937 gen(rd());
	default_random_engine generator;
	vector<Particle> newpList;
	discrete_distribution<int> d(weights.begin(), weights.end());

	for (int n = 0; n<num_particles; ++n) 
		newpList.push_back(particles[d(generator)]); //select a weighted random particle 
	
	particles = newpList;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}


