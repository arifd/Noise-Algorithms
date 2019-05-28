//////////////////////////////////////////////////////////////////////////////////////////////////
// XORSHIFT Random is the fastest random algorithm I have encountered.
// XORSHIFT RNG algorithms taken from here: https://github.com/WebDrake/xorshift
// Adapted to work as JUCE noise oscillators
//
// These algorithms originally produce unsigned 32-bit integers. Just return the variable 'output'
// if that's what you want.
// Otherwise, I divide by the the max integer value to get a float, then * 2 - 1 to ... well you know what.
//////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

class JUCERandom : public AudioIODeviceCallback
{
private:
	Random m_rnd;

public:
	void audioDeviceAboutToStart(AudioIODevice *device) override {}
	void audioDeviceStopped() override {}

	void audioDeviceIOCallback(const float **inputChannelData, int numInputChannels,
		float **outputChannelData, int numOutputChannels, int numSamples) override
	{
		for (int i = 0; i < numSamples; ++i)
		{
			float sample = generate();
			outputChannelData[0][i] = sample;
			outputChannelData[1][i] = sample;
		}
	}

	float generate()
	{
		return jmap(m_rnd.nextFloat(), 0.0f, 1.0f, -1.0f, 1.0f);
	}
};

//==============================================================================

#include <stdlib.h>
class RAND : public AudioIODeviceCallback
{
private:
public:
	RAND() {srand(101);}
	void audioDeviceAboutToStart(AudioIODevice *device) override {}
	void audioDeviceStopped() override {}

	void audioDeviceIOCallback(const float **inputChannelData, int numInputChannels,
		float **outputChannelData, int numOutputChannels, int numSamples) override
	{
		for (int i = 0; i < numSamples; ++i)
		{
			float sample = generate();
			outputChannelData[0][i] = sample;
			outputChannelData[1][i] = sample;
		}
	}

	float generate()
	{
		return (rand() / (float)RAND_MAX) * 2 - 1;
	}
};

//==============================================================================

// found here: https://stackoverflow.com/questions/26237419/faster-than-rand
class LCG : public AudioIODeviceCallback
{
private:
	unsigned int g_seed;
public:
	LCG() {g_seed = 101;}

	// Used to seed the generator.           
	void fast_srand(int seed) {g_seed = seed;}

	void audioDeviceAboutToStart(AudioIODevice *device) override {}
	void audioDeviceStopped() override {}

	void audioDeviceIOCallback(const float **inputChannelData, int numInputChannels,
		float **outputChannelData, int numOutputChannels, int numSamples) override
	{
		for (int i = 0; i < numSamples; ++i)
		{
			float sample = generate();
			outputChannelData[0][i] = sample;
			outputChannelData[1][i] = sample;
		}
	}

	float generate()
	{
		// Compute a pseudorandom integer.
		// Output value in range [0, 32767]
		g_seed = (214013*g_seed+2531011);
    	unsigned int output = (g_seed>>16)&0x7FFF;
    	return (output / 32767.0) * 2 - 1;
	}
};



//==============================================================================

// found here: https://github.com/WebDrake/xorshift
class XOR32 : public AudioIODeviceCallback
{
public:
	void audioDeviceAboutToStart(AudioIODevice *device) override {}
	void audioDeviceStopped() override {}

	void audioDeviceIOCallback(const float **inputChannelData, int numInputChannels,
		float **outputChannelData, int numOutputChannels, int numSamples) override
	{
		for (int i = 0; i < numSamples; ++i)
		{
			float sample = generate();
			outputChannelData[0][i] = sample;
			outputChannelData[1][i] = sample;
		}
	}

	float generate()
	{
		static uint32_t y = 2463534242UL;
		y^=(y<<13); y^=(y>>17);
    	uint32_t output = (y^=(y<<15));
    	return (output / 4294967296.0) * 2 - 1;
	}
};

//==============================================================================

// found here: https://jblevins.org/projects/mt
#include <cassert>
/**
 * Mersenne Twister.
 *
 * M. Matsumoto and T. Nishimura, "Mersenne Twister: A
 * 623-dimensionally equidistributed uniform pseudorandom number
 * generator", ACM Trans. on Modeling and Computer Simulation Vol. 8,
 * No. 1, January pp.3-30 (1998).
 *
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html.
 */
class MersenneTwister : public AudioIODeviceCallback
{
private:
    static const int N                    = 624;
    static const int M                    = 397;
    // constant vector a
    static const unsigned long MATRIX_A   = 0x9908b0dfUL;
    // most significant w-r bits
    static const unsigned long UPPER_MASK = 0x80000000UL;
    // least significant r bits
    static const unsigned long LOWER_MASK = 0x7fffffffUL;

    unsigned long* mt_;                  // the state vector
    int mti_;                            // mti == N+1 means mt not initialized

    unsigned long* init_key_;            // Storage for the seed vector
    int key_length_;                     // Seed vector length
    unsigned long s_;                    // Seed integer
    bool seeded_by_array_;               // Seeded by an array
    bool seeded_by_int_;                 // Seeded by an integer

public:
    MersenneTwister(void): 
	    mt_(new unsigned long[N]), mti_(N+1),
	    init_key_(NULL), key_length_(0), s_(0),
	    seeded_by_array_(false), seeded_by_int_(false)
	{
	    unsigned long init[4] = { 0x123, 0x234, 0x345, 0x456 };
	    unsigned long length = 4;
	    init_by_array(init, length);
	};

    ~MersenneTwister(void)
    {
    	assert(mt_ != NULL);
    	delete[] mt_;
    	mt_ = NULL;

    	assert(init_key_ != NULL);
    	delete[] init_key_;
    	init_key_ = NULL;
	}

    /**
	* Initializes the Mersenne Twister with a seed.
	*
 	* \param s seed
 	*/
    void init_genrand(unsigned long s)
    {
	    mt_[0]= s & 0xffffffffUL;
	    for (mti_=1; mti_<N; mti_++) {
	        mt_[mti_] = 
		    (1812433253UL * (mt_[mti_-1] ^ (mt_[mti_-1] >> 30)) + mti_); 
	        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
	        /* In the previous versions, MSBs of the seed affect   */
	        /* only MSBs of the array mt_[].                        */
	        /* 2002/01/09 modified by Makoto Matsumoto             */
	        mt_[mti_] &= 0xffffffffUL;
	        /* for >32 bit machines */
	    }
	    // Store the seed
	    s_ = s;
	    seeded_by_array_ = false;
    	seeded_by_int_ = true;
	}

	/**
	 * Seed the Mersenne Twister using an array.
	 *
	 * \param init_key an array for initializing keys
	 * \param key_length the length of \a init_key
	 */
    void init_by_array(unsigned long* init_key, int key_length)
	{
	    // Store the key array
	    int i, j, k;
	    init_genrand(19650218UL);
	    i=1; j=0;
	    k = (N>key_length ? N : key_length);
	    for (; k; k--) {
	        mt_[i] = (mt_[i] ^ ((mt_[i-1] ^ (mt_[i-1] >> 30)) * 1664525UL))
	          + init_key[j] + j; /* non linear */
	        mt_[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
	        i++; j++;
	        if (i>=N) { mt_[0] = mt_[N-1]; i=1; }
	        if (j>=key_length) j=0;
	    }
	    for (k=N-1; k; k--) {
	        mt_[i] = (mt_[i] ^ ((mt_[i-1] ^ (mt_[i-1] >> 30)) * 1566083941UL))
	          - i; /* non linear */
	        mt_[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
	        i++;
	        if (i>=N) { mt_[0] = mt_[N-1]; i=1; }
	    }

	    mt_[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 

	    // Store the seed
	    if (init_key_ != NULL) {
	        delete[] init_key_;
	    }
	    init_key_ = new unsigned long[key_length];
	    for (int k = 0; k < key_length; k++) {
	        init_key_[k] = init_key[k];
	    }
	    key_length_ = key_length;
	    seeded_by_int_ = false;
	    seeded_by_array_ = true;
	}

    /**
	* Generates a random number on [0,0xffffffff]-interval
	*
	* \return random number on [0, 0xffffffff]
 	*/
    unsigned long genrand_int32(void)
	    {
	    unsigned long y;
	    static unsigned long mag01[2]={0x0UL, MATRIX_A};
	    /* mag01[x] = x * MATRIX_A  for x=0,1 */

	    if (mti_ >= N) { /* generate N words at one time */
	        int kk;

	        if (mti_ == N+1)   /* if init_genrand() has not been called, */
	            init_genrand(5489UL); /* a default initial seed is used */

	        for (kk=0;kk<N-M;kk++) {
	            y = (mt_[kk]&UPPER_MASK)|(mt_[kk+1]&LOWER_MASK);
	            mt_[kk] = mt_[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
	        }
	        for (;kk<N-1;kk++) {
	            y = (mt_[kk]&UPPER_MASK)|(mt_[kk+1]&LOWER_MASK);
	            mt_[kk] = mt_[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
	        }
	        y = (mt_[N-1]&UPPER_MASK)|(mt_[0]&LOWER_MASK);
	        mt_[N-1] = mt_[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

	        mti_ = 0;
	    }
	  
	    y = mt_[mti_++];

	    /* Tempering */
	    y ^= (y >> 11);
	    y ^= (y << 7) & 0x9d2c5680UL;
	    y ^= (y << 15) & 0xefc60000UL;
	    y ^= (y >> 18);

	    return y;
	}

	//-------- JUCE AUDIOIODEVICECALLBACK BELOW ----- //

	void audioDeviceAboutToStart(AudioIODevice *device) override {}
	void audioDeviceStopped() override {}

	void audioDeviceIOCallback(const float **inputChannelData, int numInputChannels,
		float **outputChannelData, int numOutputChannels, int numSamples) override
	{
		for (int i = 0; i < numSamples; ++i)
		{
			float sample = generate();
			outputChannelData[0][i] = sample;
			outputChannelData[1][i] = sample;
		}
	}

	float generate()
	{
		return (genrand_int32() / 4294967295.0) * 2 - 1; 
	}
};

//==============================================================================

#include <random>
class STDMT19937 : public AudioIODeviceCallback
{
private:
	std::mt19937 mt{101};
public:
	void audioDeviceAboutToStart(AudioIODevice *device) override {}
	void audioDeviceStopped() override {}

	void audioDeviceIOCallback(const float **inputChannelData, int numInputChannels,
		float **outputChannelData, int numOutputChannels, int numSamples) override
	{
		for (int i = 0; i < numSamples; ++i)
		{
			float sample = generate();
			outputChannelData[0][i] = sample;
			outputChannelData[1][i] = sample;
		}
	}

	float generate()
	{
    	return (mt() / 4294967296.0) * 2 - 1;
	}
};

//==============================================================================

// Taken from here: https://github.com/supercollider/supercollider/blob/develop/include/plugin_interface/SC_RGen.h

//----------------------------------------------------------------------------//
// Ran088: L'Ecuyer's 1996 three-component Tausworthe generator "taus88"
//----------------------------------------------------------------------------//
//
// Returns an integer random number uniformly distributed within [0,4294967295]
//
// The period length is approximately 2^88 (which is 3*10^26).
// This generator is very fast and passes all standard statistical tests.
//
// Reference:
//   (1) P. L'Ecuyer, Maximally equidistributed combined Tausworthe generators,
//       Mathematics of Computation, 65, 203-213 (1996), see Figure 4.
//   (2) recommended in:
//       P. L'Ecuyer, Random number generation, chapter 4 of the
//       Handbook on Simulation, Ed. Jerry Banks, Wiley, 1997.
//
//----------------------------------------------------------------------------//

class Taus88 : public AudioIODeviceCallback
{
private:
	int seed = 101;
	uint32 s1, s2, s3; // random generator state
public:
	int32 Hash(int32 inKey)
	{
	    // Thomas Wang's integer hash.
	    // http://www.concentric.net/~Ttwang/tech/inthash.htm
	    // a faster hash for integers. also very good.
	    uint32 hash = (uint32)inKey;
	    hash += ~(hash << 15);
	    hash ^=   hash >> 10;
	    hash +=   hash << 3;
	    hash ^=   hash >> 6;
	    hash += ~(hash << 11);
	    hash ^=   hash >> 16;
	    return (int32)hash;
	}

	Taus88()
	{
		// humans tend to use small seeds - mess up the bits
		seed = (uint32)Hash((int)seed);

		// initialize seeds using the given seed value taking care of
		// the requirements. The constants below are arbitrary otherwise
		s1 = 1243598713U ^ seed; if (s1 <  2) s1 = 1243598713U;
		s2 = 3093459404U ^ seed; if (s2 <  8) s2 = 3093459404U;
		s3 = 1821928721U ^ seed; if (s3 < 16) s3 = 1821928721U;
	}

	uint32 trand( uint32& s1, uint32& s2, uint32& s3 )
	{
		// This function is provided for speed in inner loops where the
		// state variables are loaded into registers.
		// Thus updating the instance variables can
		// be postponed until the end of the loop.
		s1 = ((s1 &  (uint32)-2) << 12) ^ (((s1 << 13) ^  s1) >> 19);
		s2 = ((s2 &  (uint32)-8) <<  4) ^ (((s2 <<  2) ^  s2) >> 25);
		s3 = ((s3 & (uint32)-16) << 17) ^ (((s3 <<  3) ^  s3) >> 11);
		return s1 ^ s2 ^ s3;
	}

	void audioDeviceAboutToStart(AudioIODevice *device) override {}
	void audioDeviceStopped() override {}

	void audioDeviceIOCallback(const float **inputChannelData, int numInputChannels,
		float **outputChannelData, int numOutputChannels, int numSamples) override
	{
		for (int i = 0; i < numSamples; ++i)
		{
			float sample = generate();
			outputChannelData[0][i] = sample;
			outputChannelData[1][i] = sample;
		}
	}

	float generate()
	{
		return (trand(s1, s2, s3) / 4294967296.0) * 2 - 1;
	}
};
