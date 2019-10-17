template <class Type>
struct CGVectorWrapper
{
    using T = Type;
    using AggregateType = std::vector<T>;
    
    AggregateType& m_data;

    std::function<void(const AggregateType& , AggregateType& )> m_multiplyCb;
    std::function<void(AggregateType& )> m_projectCb;
    
    CGVectorWrapper(AggregateType& data)
        :m_data(data)
    {}

    AggregateType getValue()
    {
        return this->m_data;
    }
    
    CGVectorWrapper& operator = (const CGVectorWrapper& v)
    {
        m_data = v.m_data;
        return *this;
    }

    CGVectorWrapper& deepCopy(const CGVectorWrapper& v)
    {
        for(int i = 0; i < v.m_data.size(); i++) {
            m_data[i] = v.m_data[i];
        }
        return *this;
    }
    
    CGVectorWrapper& operator -= (const CGVectorWrapper& v)
    {
        for(int i = 0; i < m_data.size(); i++)
            m_data[i] -= v.m_data[i];
        return *this;
    }

    // replaces current vector with c * x + y
    void Saxpy(const T c, const CGVectorWrapper& x, const CGVectorWrapper& y)
    {
        for(int i = 0; i < m_data.size(); i++)
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
        for(int i = 0; i < m_data.size(); i++)
            result += (double) m_data[i]*(x.m_data[i]);
        return (T) result;
    }
    
    void registerMultiplyCb(std::function<void(const AggregateType& , AggregateType& )> cb) {
        if(cb)
            m_multiplyCb = cb;
    }
    
    void registerProjectCb(std::function<void(AggregateType& )> cb) {
        if(cb)
            m_projectCb = cb;
    }
    
    // replaces current vector with x
    void Multiply(const CGVectorWrapper& x)
    {
        if(m_multiplyCb) {
            m_multiplyCb(x.m_data, this->m_data);
        }
        else {
            throw std::logic_error("Multiply method not provided.");
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
