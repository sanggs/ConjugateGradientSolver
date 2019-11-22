#include <cmath>

template <class Type, int gridWidth>
struct CGVectorWrapper
{
    using T = Type;
    using AggregateType = T(&)[gridWidth];
    
    AggregateType& m_data;

    std::function<void(AggregateType& , AggregateType& )> m_multiplyCb;
    std::function<void(AggregateType& )> m_projectCb;
    
    CGVectorWrapper(AggregateType& data)
        :m_data(data)
    {}

//    AggregateType getValue()
//    {
//        return this->m_data;
//    }
    
    CGVectorWrapper& operator = (const CGVectorWrapper& v)
    {
        for(int i = 0; i < gridWidth; i++)
            m_data[i] = v.m_data[i];
        return *this;
    }

    CGVectorWrapper& deepCopy(const CGVectorWrapper& v)
    {
        for(int i = 0; i < v.m_data.size(); i++) {
            m_data[i] = v.m_data[i];
        }
        return *this;
    }
    
    CGVectorWrapper& operator += (const CGVectorWrapper& v)
    {
        for(int i = 0; i < gridWidth; i++)
            m_data[i] += v.m_data[i];
        return *this;
    }
    
    CGVectorWrapper& operator -= (const CGVectorWrapper& v)
    {
        for(int i = 0; i < gridWidth; i++)
            m_data[i] -= v.m_data[i];
        return *this;
    }

    // replaces current vector with c * x + y
    void Saxpy(const T c, const CGVectorWrapper& x, const CGVectorWrapper& y)
    {
        for(int i = 0; i < gridWidth; i++)
            m_data[i] = c * x.m_data[i] + y.m_data[i];
    }
    
    T ConvergenceNorm() {
        double result = 0.;
        for(const auto& v: m_data)
            result = std::max<double>(result, (v*v));
        return (T) std::sqrt(result);
    }
    
    T InnerProduct(const CGVectorWrapper& x) {
        double result = 0.;
        for(int i = 0; i < gridWidth; i++)
            result += (double) m_data[i]*(x.m_data[i]);
        return (T) result;
    }
    
    void registerMultiplyCb(std::function<void(AggregateType& , AggregateType& )> cb) {
        if(cb)
            m_multiplyCb = cb;
    }
    
    void registerProjectCb(std::function<void(AggregateType& )> cb) {
        if(cb)
            m_projectCb = cb;
    }
    
    // replaces current vector with x
    void MultiplyWithStencil(const CGVectorWrapper& x)
    {
        if(m_multiplyCb) {
            m_multiplyCb(x.m_data, this->m_data);
        }
        else {
            throw std::logic_error("Multiply method not provided.");
        }
    }
    
    void Multiply(const CGVectorWrapper& x) {
        for(int i = 0; i < gridWidth; i++) {
            m_data[i] = m_data[i]*(x.m_data[i]);
        }
    }
    
    void MultiplyWithScalar(const T val) {
        for(int i = 0; i < gridWidth; i++) {
            m_data[i] = m_data[i]*val;
        }
    }
    
    void ProjectToZero() {
        if (m_projectCb) {
            m_projectCb(m_data);
        }
        else {
            throw std::logic_error("Project to zero not provided.");
        }
    }
    
};
