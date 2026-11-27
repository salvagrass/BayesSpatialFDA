#pragma once

#include "SamplingUtils.hpp"
#include <functional>
#include <stdexcept>


class Recorder{
public:
    Recorder(size_t num_skip_,size_t thin_,std::string name_) 
    : num_skip(num_skip_), thin(thin_),name(name_) {}

    virtual void observe(
        const SamplerContext& ctx, 
        const ModelState& state, 
        const SamplerMemory& mem
    ) = 0;
    virtual void export_to_output(ChainOutput& out) = 0;

    virtual ~Recorder() = default;
protected:
    inline bool should_record(size_t iter) const{
        return (iter>=num_skip) && ( (iter-num_skip)% thin == 0);
    }

    inline size_t get_num_saved(size_t num_iter) const{
        return (num_iter-num_skip - 1)/thin + 1;
    }
    size_t num_skip;
    size_t thin;
    std::string name;
    size_t save_idx = 0;
};


class ScalarRecorder : public Recorder{
    using Extractor = std::function<double(const SamplerContext& ctx,const ModelState& state,const SamplerMemory& mem)>;
public:
    ScalarRecorder(size_t num_iter,size_t num_skip,size_t thin,std::string name,Extractor extractor_)
    : Recorder(num_skip,thin,name), 
    extractor(std::move(extractor_)),
    save_buf(get_num_saved(num_iter),0.0){}

    void observe(
        const SamplerContext& ctx,
        const ModelState& state,
        const SamplerMemory& mem
    ) override{
        if(this->should_record(mem.curr_iter))
            save_buf[this->save_idx++] = extractor(ctx,state,mem);
    };
    void export_to_output(ChainOutput& out) override{
        out.scalars[this->name] = std::move(this->save_buf);
    };

private:
    Extractor extractor;   
    std::vector<double> save_buf;
};


class VectorRecorder : public Recorder{
using Vector = Eigen::VectorXd;
using Extractor = std::function<Vector(const SamplerContext& ctx,const ModelState& state,const SamplerMemory& mem)>;

public:
    VectorRecorder(size_t num_iter,size_t num_skip,size_t thin,size_t vec_size,std::string name,Extractor extractor_)
    : Recorder(num_skip,thin,name), 
    extractor(std::move(extractor_)),
    save_buf(get_num_saved(num_iter),Vector::Zero(vec_size)){}

    void observe(
        const SamplerContext& ctx,
        const ModelState& state,
        const SamplerMemory& mem
    ) override {
        if(this->should_record(mem.curr_iter))
            save_buf[this->save_idx++] = extractor(ctx,state,mem);
    }
    void export_to_output(ChainOutput& out) override{
        out.vectors[this->name] = std::move(save_buf);
    }


private:
    Extractor extractor;
    std::vector<Vector> save_buf;
};

class MatrixRecorder : public Recorder{
using Matrix = Eigen::MatrixXd;
using Extractor = std::function<Matrix(const SamplerContext& ctx,const ModelState& state,const SamplerMemory& mem)>;

public:
    MatrixRecorder(size_t num_iter,size_t num_skip,size_t thin,size_t num_rows,size_t num_cols,std::string name,Extractor extractor_)
    : Recorder(num_skip,thin,name), 
    extractor(std::move(extractor_)),
    save_buf(get_num_saved(num_iter),Matrix::Zero(num_rows,num_cols)){}

    void observe(
        const SamplerContext& ctx,
        const ModelState& state,
        const SamplerMemory& mem
    ) override {
        if(this->should_record(mem.curr_iter))
            save_buf[this->save_idx++] = extractor(ctx,state,mem);
    }

    void export_to_output(ChainOutput& out) override{
        out.matrices[this->name] = std::move(save_buf);
    }

private:
    Extractor extractor;
    std::vector<Matrix> save_buf;
};