//
// Created by Steven Maio on 5/22/26.
//

#ifndef OEDLIB_LAZY_EVAL_DATA_HPP
#define OEDLIB_LAZY_EVAL_DATA_HPP

namespace OED
{
    class Lazy_Eval_Data
    {
    public:
        double gain;
        int sensor;
        int time;

        Lazy_Eval_Data(double gain, int sensor, int time)
            : gain(gain), sensor(sensor), time(time) {}

        bool operator< (const Lazy_Eval_Data &rhs) const
        {
            return gain < rhs.gain;
        }
    };
}

#endif //OEDLIB_LAZY_EVAL_DATA_HPP
