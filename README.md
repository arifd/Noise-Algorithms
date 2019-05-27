# Noise-Algorithms

When generating audible noise, the ear doesn't have a particularly high level of sophistication for detecting randomness.
Therefore we don't need the "highest quality" random number generator, but just the fastest.

In this simple console application we test a few different algorithms for noise and benchmark their speed.

It features:
- JUCE's built in Random class,
- std::rand()
- Linear congruential generator
- Xorshift
- custom implementation of Marsenne-Twister
- std::mt19937

These examples are built inside the JUCE framework, but extracting them for use anywhere is simple enough.

To compile and run this app as is, you're going to want to get and install the JUCE framework first.

https://juce.com/
