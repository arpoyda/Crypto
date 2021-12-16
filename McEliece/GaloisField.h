//
// Created by Бушев Дмитрий on 04.12.2021.
//

#ifndef MCELIECECRYPTO_GALOISFIELD_H
#define MCELIECECRYPTO_GALOISFIELD_H

#include <cinttypes>
#include <array>
#include <string>
#include <tuple>
#include <vector>
#include <optional>
#include <iostream>
#include <map>
#include <deque>
#include <cassert>
#include <random>

namespace GaloisMath {


    extern class RandomEngine {
    public:
//        RandomEngine(): rd(), gen(rd()), distrib(0u, 0xFFFFFFFF) {
//        }
        uint32_t get() {
            return std::rand();
        };
    private:
//        std::random_device rd;
//        std::mt19937 gen; //Standard mersenne_twister_engine seeded with rd()
//        std::uniform_int_distribution<> distrib;
    } randomEngine;

    uint32_t pow(uint32_t x, uint32_t degree);

    std::pair<uint32_t, uint32_t> divide_with_remainder(uint32_t num, uint32_t divisor);

    // modulo must be prime number!
    uint32_t divide_modulo(uint32_t a, uint32_t b, uint32_t modulo);
    uint32_t inverse_modulo(uint32_t a, uint32_t modulo);

    template<uint32_t p>
    class GaloisElem{
        uint32_t coef;
    public:
        GaloisElem(uint32_t num = 0u): coef(num % p){}

        explicit operator uint32_t () const{
            return coef;
        }

        GaloisElem const& operator+=(GaloisElem const& another){
            coef += another.coef;
            coef %= p;
            return *this;
        }

        GaloisElem operator+(GaloisElem const& another) const{
            GaloisElem ret{*this};
            ret += another;
            return ret;
        }

        GaloisElem const& operator-=(GaloisElem const& another){
            coef += p - another.coef;
            coef %= p;
            return *this;
        }

        GaloisElem operator-(GaloisElem const& another) const{
            GaloisElem ret{*this};
            ret -= another;
            return ret;
        }

        GaloisElem const& operator*=(GaloisElem const& another){
            coef *= another.coef;
            coef %= p;
            return *this;
        }

        GaloisElem operator*(GaloisElem const& another) const{
            GaloisElem ret{*this};
            ret *= another;
            return ret;
        }

        GaloisElem operator~() const{
            return GaloisElem{inverse_modulo(coef, p)};
        }

        GaloisElem const& operator/=(GaloisElem const& another){
            coef = divide_modulo(coef, another.coef, p);
            return *this;
        }

        GaloisElem operator/(GaloisElem const& another) const{
            GaloisElem ret{*this};
            ret /= another;
            return ret;
        }

        bool operator==(uint32_t num) const{
            return coef % p == num;
        }
        bool operator!=(uint32_t num) const{
            return coef % p != num;
        }
        bool operator==(GaloisElem const& another) const{
            return coef == another.coef;
        }
        bool operator!=(GaloisElem const& another) const{
            return coef != another.coef;
        }


        GaloisElem operator-() const{
            return GaloisElem{(p - coef) % p};
        }

    };


    template <uint32_t p>
    std::ostream& operator<<(std::ostream& stream, GaloisElem<p> elem){
        return stream << elem.operator unsigned int();
    }
    template<uint32_t m, uint32_t p>
    class FieldElemInterface{
    public:


        virtual GaloisElem<p> get(uint32_t i) const = 0;
        virtual void set(uint32_t i, GaloisElem<p> elem) = 0;

        FieldElemInterface() = default;

        bool operator==(FieldElemInterface const& another) const{
            for(int i = 0; i < m; ++i)
                if(get(i) != another.get(i))
                    return false;
            return true;
        }

        bool operator!=(FieldElemInterface const& another) const{
            return !(*this == another);
        }


        bool isNull() const{
            for(int i = m - 1; i >= 0; i--)
                if(get(i) != 0)
                    return false;
            return true;
        }

        bool neutral() const {
            for(int i = m - 1; i > 0; i--)
                if(get(i) != 0)
                    return false;
            return  get(0) == 1u;
        }


        virtual ~FieldElemInterface() = default;
    };

    template<uint32_t m, uint32_t p>
    class FieldElem: public FieldElemInterface<m, p> {
        std::array<GaloisElem<p>, m> coefs_;
    public:
        template<uint32_t m1>
        explicit FieldElem(FieldElemInterface<m1, p> const& another){
            static_assert(m1 <= m && "Cannot create field element from greater");
            for(int i = 0; i < m; i++)
                set(i, i < m1 ? another.get(i) : 0u);
        }

        template<uint32_t m1>
        FieldElem& operator=(FieldElemInterface<m1, p> const& another){
            static_assert(m1 <= m && "Cannot create field element from greater");
            for(int i = 0; i < m; i++)
                set(i, i < m1 ? another.get(i) : 0u);
        }

        GaloisElem<p> get(uint32_t i) const override{

            return coefs_.at(i);
        };
        void set(uint32_t i, GaloisElem<p> elem) override{
            coefs_.at(i) = elem;
        }


        explicit FieldElem(uint32_t num = 0u): coefs_{0u} {
            if (num > pow(p, m))
                throw std::runtime_error("Cannot create " + std::to_string(num) + " in GF(" + std::to_string(p) + "," + std::to_string(m) + ")");

            for(int i = m - 1; i >= 0; i--){
                uint32_t tmp;
                std::tie(tmp, num) = divide_with_remainder(num, pow(p, i));
                set(i, tmp);
            }
        };


        FieldElem const& operator=(uint32_t num) {
            if (num > pow(p, m))
                throw std::runtime_error("Cannot create " + std::to_string(num) + " in GF(" + std::to_string(p) + "," + std::to_string(m) + ")");

            for(int i = m - 1; i >= 0; i--){
                uint32_t tmp;
                std::tie(tmp, num) = divide_with_remainder(num, pow(p, i));
                set(i, tmp);
            }

            return *this;
        }


        FieldElem const &operator+=(FieldElem const &another){
            for(int i = 0; i < m; ++i){
                coefs_.at(i) += another.coefs_.at(i);
            }

            return *this;
        }

        FieldElem operator+(FieldElem const &another) const{
            FieldElem ret{*this};
            ret += another;
            return ret;
        }

        FieldElem const &operator-=(FieldElem const &another){
            for(int i = 0; i < m; ++i){
                coefs_.at(i) -= another.coefs_.at(i);
            }

            return *this;
        }
        FieldElem operator-(FieldElem const &another) const{
            FieldElem ret{*this};
            ret -= another;
            return ret;
        }

        FieldElem operator-() const{
            FieldElem ret;
            for(int i = 0; i < m; ++i){
                ret.coefs_.at(i) = - coefs_.at(i);
            }

            return ret;
        }

        uint32_t get() const{
            uint32_t ret = 0u;
            for(int i = 0; i < m; i++)
                ret += static_cast<uint32_t>(get(i)) * pow(p, i);

            return ret;
        }

        void print() const {
            std::cout << "<" << p << "," << m << ">: ";
            bool printedFlag = false;
            for(int i = m - 1; i >= 0; --i ) {
                if(get(i) == 0)
                    continue;
                if(i != m - 1 && printedFlag)
                    std::cout << " + ";
                if(get(i) != 1 || i == 0)
                    std::cout << get(i);
                if(i != 0)
                    std::cout << "x^" << i;

                printedFlag = true;
            }
            std::cout << " = " << get() << std::endl;
        }

    };

    constexpr uint32_t pow_const(uint32_t x, uint32_t degree){
        if(degree == 0)
            return 1u;
        if(degree == 1)
            return x;
        uint32_t lhs = degree / 2u;
        uint32_t rhs = degree - lhs;

        return pow_const(x, lhs) * pow_const(x, rhs);
    }

    constexpr uint32_t modulo_digit(uint32_t num, uint32_t modulo, uint32_t digit)
    {
        if(digit == 0)
            return num % modulo;
        for(int i = 0; i < digit; i++){
            num -= num % pow_const(modulo, i + 1);
        }
        return num / pow_const(modulo, digit) % modulo;
    }

    template<uint32_t num, uint32_t N, uint32_t... Rest>
    struct GaloisIrreducible_impl {
        static constexpr auto& value = GaloisIrreducible_impl<num,  N - 1, modulo_digit(num, 2u, N), Rest...>::value;
    };

    template<uint32_t num, uint32_t... Rest>
    struct GaloisIrreducible_impl<num, 0, Rest...> {
        static constexpr uint32_t value[] = { modulo_digit(num, 2u, 0), Rest... };
    };


    template<uint32_t num,uint32_t... Rest>
    constexpr uint32_t GaloisIrreducible_impl<num, 0, Rest...>::value[];

    template<uint32_t m, uint32_t num>
    class GaloisIrreducible {
        static_assert(m >= 0, "N must be at least 0");
    public:
        static constexpr auto& value = GaloisIrreducible_impl<num, m>::value;

        static const struct Accessor: public FieldElemInterface<m + 1, 2>{
            GaloisElem<2> get(uint32_t i) const override{
                uint32_t b;
                if(i > m)
                    throw std::runtime_error("Out of bounds (size = " + std::to_string(m) + ", i = " + std::to_string(i) + ")");
                return value[i];
            };
            void set(uint32_t i, GaloisElem<2> elem) override{
            };
        } accessor;

        GaloisIrreducible() = delete;
        GaloisIrreducible(const GaloisIrreducible&) = delete;
        GaloisIrreducible(GaloisIrreducible&&) = delete;
    };
    template<uint32_t m, uint32_t num>
    typename GaloisIrreducible<m, num>::Accessor const GaloisIrreducible<m, num>::accessor;

    template<typename U, uint32_t m1, uint32_t m2, uint32_t p>
    U divide(FieldElemInterface<m1, p>& divisor, FieldElemInterface<m2, p> const& divider){
        U ret{0};
        std::optional<uint32_t> highestDegree1;
        std::optional<uint32_t> highestDegree2;
        for(int i = m1 - 1; i >= 0; i--)
            if(divisor.get(i) != 0u) {
                highestDegree1 = i;
                break;
            }

        for(int i = m2 - 1; i >= 0; i--)
            if(divider.get(i) != 0u) {
                highestDegree2 = i;
                break;
            }

        if(!highestDegree2)
            throw std::runtime_error("Division by null polynomial");

        if(!highestDegree1){
            return ret;
        }

        auto cachedHighDegree2 = highestDegree2.value();

        for(int i = highestDegree1.value(); i >= 0; i--){
            if(!highestDegree1)
                break;

            auto cachedHighDegree = highestDegree1.value();

            int  degreeDiff = cachedHighDegree - cachedHighDegree2;
            if(degreeDiff < 0)
                break;
            auto sub = divisor.get(cachedHighDegree) / divider.get(cachedHighDegree2);
            ret.set(degreeDiff, sub);

            for(int j = cachedHighDegree2; j >= 0; j--){
                auto tmp = divisor.get(j + degreeDiff);
                tmp -= (sub * divider.get(j));
                divisor.set(j + degreeDiff, tmp);
            }
            highestDegree1.reset();
            for(int j = cachedHighDegree; j >= 0; --j){
                if(divisor.get(j) != 0u) {
                    highestDegree1.emplace(j);
                    break;
                }
            }
        }

        return ret;
    }
    template<uint32_t m, uint32_t m1, uint32_t p>
    FieldElem<m, p> mod(FieldElemInterface<m, p> const& lhs, FieldElemInterface<m1, p> const& rhs){
        FieldElem<m, p> ret{lhs};
        divide<FieldElem<m1, p>>(ret, rhs);
        return ret;
    }

    template<uint32_t m, typename Irreducible = GaloisIrreducible<m, pow_const(2u, m) + 1> /* TODO: enable only if poly degree is not over m */>
    class GaloisNumber : public FieldElem< m, 2>{
    public:
        using FieldElem<m, 2>::get;
        using FieldElem<m, 2>::set;
        using FieldElemInterface<m, 2>::neutral;
        using FieldElemInterface<m, 2>::isNull;
        GaloisNumber(uint32_t num = 0u): FieldElem<m, 2>(num) {
        };

        template<uint32_t m1>
        GaloisNumber(FieldElemInterface<m1, 2> const& another): FieldElem<m, 2>(another){}

        GaloisNumber const& operator=(uint32_t num){
            FieldElem<m, 2>::operator=(num);
            return *this;
        }

        bool operator==(GaloisNumber const& another) const{
            for(uint32_t i = 0; i < m; ++i){
                if(get(i) != another.get(i))
                    return false;
            }
            return true;
        }

        bool operator!=(GaloisNumber const& another) const{
            return !(*this == another);
        }


        GaloisNumber const &operator*=(GaloisNumber const &another){

            FieldElem<2u * m, 2> multiplied{0u};
            for(int i = 0; i < m; i++)
                for(int j = 0; j < m; j++)
                    multiplied.set(i + j, (multiplied.get(i + j) + get(i) * another.get(j)));

            divide<FieldElem<2u * m, 2>>(multiplied, Irreducible::accessor);


            for(int i = 0; i < m; i++)
                set(i, multiplied.get(i));
            return *this;
        }
        GaloisNumber operator*(GaloisNumber const &another) const{
            GaloisNumber ret = *this;
            ret *= another;
            return ret;
        }

        GaloisNumber operator^(uint32_t degree) const{
            GaloisNumber ret{1u};
            for(auto i = 0u; i < degree; ++i)
                ret *= *this;
            return ret;
        }

        GaloisNumber const &operator/=(GaloisNumber const &another){
            *this*=~another;
            return *this;
        }
        GaloisNumber operator/(GaloisNumber const &another) const{
            GaloisNumber ret = *this;
            ret /= another;
            return ret;
        }


        // inverse
        GaloisNumber operator~() const{
            if(neutral())
                return *this;

            std::vector<std::pair<GaloisNumber<m + 1>, GaloisNumber<m + 1>>> euclid_stack;

            // euclid algorithm
            GaloisNumber<m + 1> greater{Irreducible::accessor};
            GaloisNumber<m + 1> smaller = *this;

            euclid_stack.emplace_back(0u, greater);
            euclid_stack.emplace_back(0u, smaller);

            while (!smaller.neutral()){
                auto division = divide<GaloisNumber<m + 1>>(greater, smaller);
                euclid_stack.emplace_back(division, greater);
                auto tmp = greater;
                greater = smaller;
                smaller = tmp;
            }

            GaloisNumber<m + 1, Irreducible> q, ac, p_, bc;

            q = greater;
            p_ = smaller;
            ac = 1;
            bc = -euclid_stack.back().first;


            for(int i = euclid_stack.size() - 1; i >= 3; --i){

                auto tmp_ac = ac;
                ac = bc;
                bc = tmp_ac - euclid_stack.at(i - 1).first * bc;
                p_ = q;
                q = euclid_stack.at(i - 3).second;

            }

            GaloisNumber ret;

            for(int i = 0; i < m; i++)
                ret.set(i, bc.get(i));
            return ret;
        }

        static GaloisNumber generateIrreducible(uint32_t degree){
            GaloisNumber ret{0u};
            ret.set(degree, 1u);
            ret.set(0, 1u);

            uint32_t degreePow = pow(2u, degree);
            uint32_t totalElements = pow(2u, m);
            for(uint32_t i = 1; i < degreePow; i++){
                auto tmp = ret + GaloisNumber{i};
                bool primeCheck = true;
                for(uint32_t j = 2; j < totalElements; ++j){
                    if(j == i + degreePow + 1u)
                        continue;
                    auto tmp2 = tmp;
                    divide<GaloisNumber>(tmp2, GaloisNumber{j});
                    if(tmp2.isNull()){
                        primeCheck = false;
                        break;
                    }
                }
                if(primeCheck)
                    return tmp;
            }
            throw std::runtime_error("Could not generate irreducible of degree = " +
                    std::to_string(degree) + " in galois field <" +
                    std::to_string(2u) + "," + std::to_string(m) + ">");

            return ret;
        }
    };

    template<uint32_t m, typename U>
    std::ostream& operator<<(std::ostream& stream, GaloisNumber<m, U> num){
        num.print();
        return stream;
    }


    template<uint32_t m, uint32_t gf>
    class PolynomialOverExtendedField{
    public:
        using Member = GaloisNumber<m, GaloisIrreducible<m, gf>>;
    private:
        std::vector<Member> coefs_;
        void fit(){
            if(!coefs_.empty()) {
                uint32_t leadingCoef = coefs_.size() - 1;
                for (int i = leadingCoef; i >= 0; --i) {
                    if (!coefs_.at(i).isNull()) {
                        coefs_.resize(i + 1);
                        return;
                    }
                }
            }
            coefs_.clear();
            coefs_.emplace_back(0u);
        }

    public:
        void print() const{
            for(uint32_t i = 0; i < coefs_.size(); ++i){
                std::cout << "X^" << i << ": ";
                coefs_.at(i).print();
            }
        }
        explicit PolynomialOverExtendedField(std::vector<uint32_t> const& coefs = {}, bool b = true){
            for(auto const& c: coefs){
                coefs_.push_back(Member{2}^c);
            }
            fit();
            //std::cout << "Goppa polynomial:" << std::endl;
            //print();
        };

        explicit PolynomialOverExtendedField(std::vector<Member> coefs, int b = 0): coefs_(coefs){
            fit();
            //std::cout << "Goppa polynomial:" << std::endl;
            //print();
        };

        void set(uint32_t i, Member member) {
            if(i >= coefs_.size()){
                coefs_.resize(i + 1);
            }

            coefs_.at(i) = member;
            fit();
        }

        uint32_t degree() const{
            return coefs_.size() - 1;
        }

        void clear(){
            coefs_.clear();
            fit();
        }

        bool isNull() const{
            return degree() == 0 && coefs_.at(0).isNull();
        }

        bool isNeutral() const {
            return degree() == 0 && (coefs_.at(0) == 1u);
        }



        PolynomialOverExtendedField const& operator+=(PolynomialOverExtendedField const& another){
            if(another.degree() > degree()){
                coefs_.resize(another.coefs_.size());
            }

            for(uint32_t i = 0u; i < another.coefs_.size(); ++i)
                coefs_.at(i) += another.coefs_.at(i);

            fit();
            return *this;
        }

        PolynomialOverExtendedField operator+(PolynomialOverExtendedField const& another) const{
            PolynomialOverExtendedField ret{*this};
            ret += another;
            return ret;
        }

        PolynomialOverExtendedField const& operator-=(PolynomialOverExtendedField const& another){
            if(another.degree() > degree()){
                coefs_.resize(another.coefs_.size());
            }
            for(uint32_t i = 0u; i < another.coefs_.size(); ++i)
                coefs_.at(i) -= another.coefs_.at(i);

            fit();
            return *this;
        }

        PolynomialOverExtendedField operator-(PolynomialOverExtendedField const& another) const{
            PolynomialOverExtendedField ret{*this};
            ret -= another;
            return ret;
        }

        PolynomialOverExtendedField operator-() const{
            PolynomialOverExtendedField ret;
            ret -= *this;
            return ret;
        }

        bool operator==(PolynomialOverExtendedField const& another) const{
            if(degree() != another.degree())
                return false;
            for(uint32_t i = 0; i < coefs_.size(); ++i)
                if(coefs_.at(i) != another.coefs_.at(i))
                    return false;

            return true;
        }

        bool operator!=(PolynomialOverExtendedField const& another) const{
            return !(*this == another);
        }

        PolynomialOverExtendedField divideWithRemainder(PolynomialOverExtendedField const& another){
            PolynomialOverExtendedField ret{};
            try{

            if(another.isNull())
                throw std::runtime_error("Division by null polynomial");

            uint32_t highestDegree1 = degree();
            uint32_t highestDegree2 = another.degree();


            for(int i = highestDegree1; i >= 0; i--){
                int  degreeDiff = highestDegree1 - highestDegree2;

                if(degreeDiff < 0)
                    break;

                auto sub = coefs_.at(highestDegree1) /  another.coefs_.at(highestDegree2);
                ret.set(degreeDiff, sub);

                for(int j = highestDegree2; j >= 0; j--){
                    auto tmp = coefs_.at(j + degreeDiff);
                    tmp -= (sub * another.coefs_.at(j));
                    coefs_.at(j + degreeDiff) = tmp;
                }
                fit();
                highestDegree1 = degree();
            }

            ret.fit();

            } catch(std::out_of_range const& e){
                std::cout << e.what() << std::endl;
                exit(0);
            }
            return ret;
        }

        PolynomialOverExtendedField multiply(PolynomialOverExtendedField const& another, PolynomialOverExtendedField const& modulo) const{
            PolynomialOverExtendedField multiplied;
            try {
                uint32_t deg = degree();
                uint32_t deg2 = another.degree();

                multiplied.coefs_.resize(coefs_.size() + another.coefs_.size());
                for (uint32_t i = 0; i <= deg; ++i)
                    for (uint32_t j = 0; j <= deg2; ++j)
                        multiplied.coefs_.at(i + j) += coefs_.at(i) * another.coefs_.at(j);

                multiplied.fit();
                multiplied.divideWithRemainder(modulo);

                multiplied.fit();
            } catch(std::out_of_range const& e){
                std::cout << e.what() << std::endl;
                exit(0);
            }

            return multiplied;
        }

        PolynomialOverExtendedField pow(uint32_t x, PolynomialOverExtendedField const& modulo) const{
            if(x == 0){
                PolynomialOverExtendedField ret;
                ret.coefs_.at(0) = 1u;
                return ret;
            }
            if(x == 1){
                return *this;
            }

            return pow(x / 2, modulo).multiply(pow(x - x / 2, modulo), modulo);

        }

        PolynomialOverExtendedField inverse(PolynomialOverExtendedField const& modulo) const{
            auto tmp = *this;
            PolynomialOverExtendedField accum{{1u}, 1};
            for(uint32_t i = 0; i < m * modulo.degree(); ++i){
                if(i != 0)
                    accum = accum.multiply(tmp, modulo);
                tmp = tmp.multiply(tmp, modulo);
            }

            return accum; //pow(::GaloisMath::pow(2u, m * modulo.degree()) - 2, modulo);
        }

        PolynomialOverExtendedField divide(PolynomialOverExtendedField const& another, PolynomialOverExtendedField const& modulo) const{
            return multiply(another.inverse(modulo), modulo);
        }

        Member operator()(Member const& member) const{
            auto ret = coefs_.at(0);
            auto deg = degree();
            for(uint32_t i = 1; i <= deg; ++i){
                ret += coefs_.at(i) * (member^i);
            }

            return ret;
        }

        Member g(uint32_t i) const{
            if(i >= coefs_.size())
                return Member {0u};
            return coefs_.at(i);
        }

        bool isIrreducable() const{
            return true;

            auto deg = degree();

            for(uint32_t i = 1; i < deg; ++i){
                uint32_t j = deg - i;
                PolynomialOverExtendedField a;
                PolynomialOverExtendedField b;
                a.set(i, 1u);
                b.set(j, 1u);

            }
            for(uint32_t i = 0; i < pow_const(2, m); ++i)
                if(operator()(Member{i}).isNull()) {
                    std::cout << "Root:" << std::endl;
                    Member{i}.print();
                    return false;
                }
            return true;
        }
    };

    constexpr uint32_t min(uint32_t first, uint32_t second){
        return first < second ? first : second;
    }
    template <unsigned dim1, unsigned dim2, typename T = float>
    class Matrix{
        std::vector<T> elems_{0};
    public:

        explicit Matrix(T elem = 0){
            elems_.resize(dim1 * dim2);
            if(dim1 == dim2)
                for(uint32_t i = 0; i < dim2; ++i)
                    at(i, i) = elem;
        }

        T& at(int i, int j){
            return elems_.at(j * dim1 + i);
        }

        T const& at(int i, int j) const{
            return elems_.at(j * dim1 + i);
        }
        template<unsigned adim1, unsigned adim2, typename U>
        std::enable_if_t<adim2 == dim1, Matrix<adim1, dim2, T>> operator*(Matrix<adim1, adim2, U> const& another) const{
            Matrix<adim1, dim2, T> ret;
            for(int i = 0; i < adim1; i++)
                for(int j = 0; j < dim2; j++)
                    for(int k = 0; k < dim1; k++)
                        ret.at(i, j) += at(k, j) * another.at(i, k);
            return ret;
        }

        bool isNull() const{
            for(uint32_t i = 0; i < elems_.size(); ++i)
                if(elems_.at(i) != 0)
                    return false;

            return true;
        }

        bool operator==(Matrix<dim1, dim2, T> const& another) const{
            for(int i = 0; i < dim1; i++)
                for(int j = 0; j < dim2; j++)
                    if(at(i, j) != another.at(i, j))
                        return false;
            return true;
        }

        bool operator!=(Matrix<dim1, dim2, T> const& another) const{
            return !(*this == another);
        }
        Matrix operator-() const{
            Matrix ret{*this};
            for(auto& elem: ret.elems_)
                elem = -elem;
            return ret;
        }
        template<typename U>
        Matrix const& operator-=(Matrix<dim1, dim2, U> const& another){
            for(uint32_t i = 0; i < dim2; ++i)
                for(uint32_t j = 0; j < dim1; ++j)
                    at(j, i) -= another.at(j, i);

            return *this;
        }
        template<typename U>
        Matrix operator-(Matrix<dim1, dim2, U> const& another) const{
            Matrix ret{*this};
            ret -= another;

            return ret;
        }
        template<typename U>
        Matrix const& operator+=(Matrix<dim1, dim2, U> const& another){
            for(uint32_t i = 0; i < dim2; ++i)
                for(uint32_t j = 0; j < dim1; ++j)
                    at(j, i) += another.at(j, i);

            return *this;
        }
        template<typename U>
        Matrix operator+(Matrix<dim1, dim2, U> const& another) const{
            Matrix ret{*this};
            ret += another;

            return ret;
        }
        void print() const{
            for(int i = 0; i < dim2; i++){
                for(int j = 0; j < dim1; j++)
                    std::cout << at(j,i) << ", ";

                std::cout << std::endl;
            }
        }

        template<uint32_t sub1, uint32_t sub2>
        std::enable_if_t<sub1 <= dim1 && sub2 <= dim2, Matrix<sub1, sub2, T>> submatrix(uint32_t offsetX, uint32_t offsetY) const{
            Matrix<sub1, sub2, T> ret;
            for(uint32_t i = offsetY; i < offsetY + sub2; ++i)
                for(uint32_t j = offsetX; j < offsetX + sub1; ++j)
                    ret.at(j - offsetX, i - offsetY) = at(j, i);

            return ret;
        }

        template<uint32_t sub1, uint32_t sub2>
        std::enable_if_t<sub1 <= dim1 && sub2 <= dim2, Matrix const&> copySubmatrix(Matrix<sub1, sub2, T> matrix, uint32_t offsetX, uint32_t offsetY){
            for(uint32_t i = offsetY; i < offsetY + sub2; ++i)
                for(uint32_t j = offsetX; j < offsetX + sub1; ++j)
                    at(j, i) = matrix.at(j - offsetX, i - offsetY);

            return *this;
        }
        void swapRows(uint32_t i, uint32_t j){
            if(i == j)
                return;
            for(int k = 0; k < dim1; ++k){
                auto tmp = at(k, i);
                at(k, i) = at(k, j);
                at(k, j) = tmp;
            }
        }
        template<typename U>
        void multiplyRow(U elem, uint32_t j){
            for(int k = 0; k < dim1; ++k) {
                at(k, j) *= elem;
            }
        }
        template<typename U>
        void divideRow(U elem, uint32_t j){
            for(int k = 0; k < dim1; ++k) {
                at(k, j) /= elem;
            }
        }
        template<typename U>
        void subtractRows(uint32_t i, uint32_t j, U mult){
            for(int k = 0; k < dim1; ++k) {
                at(k, i) -= at(k, j) * mult;
            }
        }

        template<typename U>
        void addRows(uint32_t i, uint32_t j, U mult){
            for(int k = 0; k < dim1; ++k) {
                at(k, i) += at(k, j) * mult;
            }
        }

        void swapColumns(uint32_t i, uint32_t j){
            for(uint32_t k = 0; k < dim2; ++k){
                auto tmp = at(i, k);
                at(i, k) = at(j, k);
                at(j, k) = tmp;
            }
        }

        Matrix<dim2, dim1, T> transpose() const{
            Matrix<dim2, dim1, T> ret;

            for(uint32_t i = 0; i < dim2; i++)
                for(uint32_t j = 0; j < dim1; j++)
                    ret.at(i, j) = at(j, i);
            return ret;
        }

        Matrix reduce() const{
            Matrix mat = *this;
            for(uint32_t i = 0; i < min(dim1, dim2); ++i){
                bool foundNonNull = false;
                for(uint32_t k = i; k < dim1; k++) {
                    for (uint32_t j = i; j < dim2; ++j)
                        if (mat.at(k, j) != 0u) {
                            if(i != k) {
                                mat.swapColumns(i, k);
                            }
                            if(i != j) {
                                mat.swapRows(i, j);
                            }
                            foundNonNull = true;
                            break;
                        }
                    if(foundNonNull)
                        break;
                }
                if(!foundNonNull)
                    return mat;

                auto div = mat.at(i, i);
                mat.divideRow(div, i);

                for(uint32_t j = 0; j < dim2; ++j){
                    if(i == j)
                        continue;
                    auto sub = mat.at(i, j);
                    mat.subtractRows(j, i, sub);
                }
            }
            return mat;
        }



        template<uint32_t adim2>
        Matrix<dim1, adim2, T> nullspace() const{
            if(dim1 < dim2)
                throw std::runtime_error("Nullspace allowed only on cols > rows matrices");
            if(dim1 - rank() != adim2)
                throw std::runtime_error("Nullspace matrix must have rows exactly as dimensions of nullspace");
            auto tmp = *this;

            std::array<uint32_t, dim1> rowsIds;
            for(uint32_t i = 0; i < dim1; ++i)
                rowsIds.at(i) = i;

            for(int i = 0; i < dim2; ++i){
                bool foundNonNull = false;
                for(uint32_t j = i; j < dim2; ++j){
                    if(tmp.at(i, j) != 0u){
                        tmp.swapRows(i, j);
                        foundNonNull = true;
                        break;
                    }
                }
                if(!foundNonNull) {
                    bool foundNotNullCol = false;
                    for(uint32_t j = i + 1; j < dim1; ++j){
                        for(uint32_t k = i; k < dim2; ++k){
                            if(tmp.at(j, k) != 0){
                                tmp.swapColumns(i, j);
                                auto idTmp = rowsIds.at(i);
                                rowsIds.at(i) = rowsIds.at(j);
                                rowsIds.at(j) = idTmp;
                                foundNotNullCol = true;
                                break;
                            }
                        }
                    }
                    if(foundNotNullCol){
                        --i;
                        continue;
                    } else{
                        break;
                    }
                }
                tmp.divideRow(tmp.at(i, i), i);
                for(uint32_t j = 0; j < dim2; ++j){
                    if(j == i)
                        continue;
                    tmp.subtractRows(j, i, tmp.at(i, j) / tmp.at(i, i));
                }
            }

            Matrix<adim2, dim2, T> almostNullspace = -tmp.template submatrix<adim2, dim2>(dim1 - adim2, 0);

            //almostNullspace.print();
            Matrix<adim2, dim1, T> nullspaceMat;
            for(uint32_t i = 0; i < adim2; ++i)
                nullspaceMat.at(i, i + dim1 - adim2) = 1;

            nullspaceMat.template copySubmatrix(almostNullspace, 0, 0);


            for(uint32_t i = 0; i < dim1; ++i)
                for(uint32_t j = i; j < dim1; ++j){
                    if(rowsIds.at(j) == i){
                        if(j == i)
                            continue;
                        auto tmpId = rowsIds.at(i);
                        rowsIds.at(i) = rowsIds.at(j);
                        rowsIds.at(j) = tmpId;
                        nullspaceMat.swapRows(i, j);
                    }
                }


            auto nullmat = nullspaceMat.transpose() * transpose();

            return nullspaceMat.transpose();

        }

        //template<std::enable_if<dim1 >= dim2, bool> = true>
        static uint32_t rankImpl(Matrix matrix){
            auto& tmp = matrix;
            uint32_t ret = 0u;

            for(int i = 0; i < dim2; ++i){
                bool foundNonNull = false;
                for(uint32_t j = i; j < dim2; ++j){
                    if(tmp.at(i, j) != 0u){
                        tmp.swapRows(i, j);
                        foundNonNull = true;
                        break;
                    }
                }
                if(!foundNonNull) {
                    bool foundNotNullCol = false;
                    for(uint32_t j = i + 1; j < dim1; ++j){
                        for(uint32_t k = i; k < dim2; ++k){
                            if(tmp.at(j, k) != 0){
                                tmp.swapColumns(i, j);
                                foundNotNullCol = true;
                                break;
                            }
                        }
                    }
                    if(foundNotNullCol){
                        --i;
                        continue;
                    } else{
                        break;
                    }

                }
                ret++;
                tmp.divideRow(tmp.at(i, i), i);
                for(uint32_t j = 0; j < dim2; ++j){
                    if(j == i)
                        continue;

                    tmp.subtractRows(j, i, tmp.at(i, j) / tmp.at(i, i));
                }
            }

            return ret;

        }
        uint32_t rank() const{
            auto tmp = *this;
            uint32_t dim1_, dim2_;
            if(dim1 > dim2){
                return rankImpl(*this);
            } else{
                return Matrix<dim2, dim1, T>::rankImpl(transpose());
            }
         }
    };


    template<typename T,uint32_t rank>
    Matrix<rank, rank, T> randomPermutation(){
        Matrix<rank, rank, T> ret;
        std::deque<uint32_t> ids;
        for(uint32_t i = 0; i < rank; ++i){
            ids.emplace_back(i);
        }

        for(uint32_t i = 0; i < rank; ++i){
            uint32_t randPick = randomEngine.get() % (rank - i);
            ret.at(i, ids.at(randPick)) = 1u;
            ids.erase(ids.begin() + randPick);
        }

        return ret;
    }

    template<typename T,uint32_t rank>
    Matrix<rank, rank, T> inverse(Matrix<rank, rank, T> const& matrix){
        Matrix<rank, rank, T> mat = matrix;
        Matrix<rank, rank, T> ret{1u};
        for(uint32_t i = 0; i < rank; ++i){
            bool foundNonNull = false;
            for(uint32_t k = i; k < rank; k++) {
                for (uint32_t j = i; j < rank; ++j)
                    if (mat.at(k, j) != 0u) {
                        if(i != k) {
                            mat.swapColumns(i, k);
                            ret.swapColumns(i, k);
                        }
                        if(i != j) {
                            mat.swapRows(i, j);
                            ret.swapRows(i, j);
                        }
                        foundNonNull = true;
                        break;
                    }
                if(foundNonNull)
                    break;
            }
            if(!foundNonNull)
                throw std::runtime_error("Matrix is not invertible");

            auto div = mat.at(i, i);
            mat.divideRow(div, i);
            ret.divideRow(div, i);

           for(uint32_t j = 0; j < rank; ++j){
               if(i == j)
                   continue;
               auto sub = mat.at(i, j);
               mat.subtractRows(j, i, sub);
               ret.subtractRows(j, i, sub);
           }
        }

        return ret;
    }
    template<typename T,uint32_t rank>
    Matrix<rank, rank, T> randomInvertible(){
        Matrix<rank, rank, T> ret{1};

        for(uint32_t i = 0; i < rank * rank; ++i){
            uint32_t randPick1 = randomEngine.get() % rank;
            uint32_t randPick2 = randomEngine.get() % rank;
            if(randPick2 == randPick1)
                continue;

            T randomPickMult = randomEngine.get();
            while(randomPickMult == 0)
                randomPickMult = randomEngine.get();
            ret.addRows(randPick1, randPick2, randomPickMult);
        }

        for(uint32_t i = 0; i < rank; ++i){
            uint32_t randPick1 = randomEngine.get() % rank;
            uint32_t randPick2 = randomEngine.get() % rank;
            ret.swapRows(randPick1, randPick2);
        }

        return ret;
    }

    template<uint32_t m, uint32_t t, uint32_t gf>
    class LinearGoppaCodeGen{
        static constexpr const uint32_t n = pow_const(2u, m);
        static constexpr const uint32_t k = n - m * t;
    public:
        using GPoly = PolynomialOverExtendedField<m, gf>;
    private:
        GPoly goppaPoly;
        Matrix<n, k, GaloisElem<2u>> g_mat_;
        Matrix<k, k, GaloisElem<2u>> S_;
        Matrix<n, n, GaloisElem<2u>> P_;
        Matrix<n, k, GaloisElem<2u>> G_public;
        std::array<typename GPoly::Member, n> alphas_;
        Matrix<n, t * m, GaloisElem<2u>> H_;
        Matrix<n, 2 * t * m, GaloisElem<2u>> H2_;

    public:
        using EncryptedMessage = Matrix<n, 1, GaloisElem<2u>>;
        using Message = Matrix<k, 1, GaloisElem<2u>>;
        using GoppaCoefs = std::vector<uint32_t>;
        using Key = std::pair<Matrix<n, k, GaloisElem<2u>>, uint32_t>;
        LinearGoppaCodeGen(GoppaCoefs const& goppa): goppaPoly(goppa){
            if(t != goppaPoly.degree()){
                throw std::runtime_error("Goppa polynom degree(which is " + std::to_string(goppaPoly.degree()) + ") must match t parameter(which is " + std::to_string(t) + ")");
            }
            uint32_t i = 0;

            if(!goppaPoly.isIrreducable())
                throw std::runtime_error("Given goppa polynomial isn't irreducible");

            alphas_.at(0) = 0u;

            for(uint32_t i = 1; i < n; ++i){
                alphas_.at(i) = typename GPoly::Member{2}^(i - 1) ;
            }



            std::sort(alphas_.begin(), alphas_.end(), [](typename GPoly::Member const& lhs, typename GPoly::Member const& rhs){ return lhs.get() < rhs.get();});

            for(uint32_t i = 1; i < n; ++i){
                if(alphas_.at(i) == alphas_.at(i - 1)){
                    throw std::runtime_error("GF<2," + std::to_string(m) + "> produced by not irreducible.");
                }
            }

            alphas_.at(0) = 0u;

            for(uint32_t i = 1; i < n; ++i){
                alphas_.at(i) = typename GPoly::Member{2}^(i - 1) ;
            }


            Matrix<t, t, typename GPoly::Member> X;

            for(i = 0; i < t; ++i)
                for(uint32_t j = i; j < t; ++j)
                    X.at(j, i) = -goppaPoly.g(t - (j - i));



            Matrix<n, t, typename GPoly::Member> Y;

            for(i = 0; i < t; ++i)
                for(uint32_t j = 0; j < n; ++j)
                    Y.at(j, i) = (alphas_.at(j)^(t - 1 - i));


            Matrix<n, n, typename GPoly::Member> Z;

            for(i = 0; i < n; ++i)
                Z.at(i, i) = ~(goppaPoly(alphas_.at(i)));



            auto HGalois = X * Y * Z;

            for(i = 0; i < t; ++i)
                for(uint32_t j = 0; j < n; ++j)
                    for(uint32_t k_ = 0; k_ < m; k_++)
                        H_.at(j, i * m + k_) = HGalois.at(j, i).get(k_);



            g_mat_ = H_.template nullspace<k>();


            assert(( g_mat_ * H_.transpose()).isNull() && "GHt must be null");


            S_ = randomInvertible<GaloisElem<2>, k>();
            P_ = randomPermutation<GaloisElem<2>, n>();

            auto eMat = Matrix<k,k, GaloisElem<2>>{1};
            assert(S_ * inverse(S_) == inverse(S_) * S_ && S_ * inverse(S_) == eMat);


            G_public = S_ * g_mat_ * P_;

            PolynomialOverExtendedField<m, gf> bigModule;
            bigModule.set(t * 3, 1u);
            auto goppaQuad = goppaPoly.pow(2, bigModule);


            Matrix<2 * t, 2 * t, typename GPoly::Member> equationMat;
            Matrix<1, 2 * t, typename GPoly::Member> freeMembers;



            Matrix<2 * t, 2 * t, typename GPoly::Member> X2;

            for(uint32_t i = 0; i < 2 * t; ++i)
                for(uint32_t j = i; j < 2 * t; ++j)
                    X2.at(j, i) = -goppaQuad.g(2 * t - (j - i));



            Matrix<n, 2 * t, typename GPoly::Member> Y2;

            for(uint32_t i = 0; i < 2 * t; ++i)
                for(uint32_t j = 0; j < n; ++j)
                    Y2.at(j, i) = (alphas_.at(j)^(2 * t - 1 - i));


            Matrix<n, n, typename GPoly::Member> Z2;

            for(uint32_t i = 0; i < n; ++i)
                Z2.at(i, i) = ~(goppaQuad(alphas_.at(i)));



            auto HGalois2 = X2 * Y2 * Z2;



            for(uint32_t i = 0; i < 2 * t; ++i)
                for(uint32_t j = 0; j < n; ++j)
                    for(uint32_t k_ = 0; k_ < m; k_++)
                        H2_.at(j, i * m + k_) = HGalois2.at(j, i).get(k_);

        }

        Matrix<n, m, GaloisElem<2>>const& G() {return g_mat_;};

        Matrix<n, m, GaloisElem<2>>const& GPub() {return G_public;};

        Key getPublicKey() const { return {G_public, t};}

        Message decrypt(EncryptedMessage const& msg){
            using Polynomial = PolynomialOverExtendedField<m, gf>;
            EncryptedMessage y = msg * P_.transpose();



            Polynomial pSyndrome;
            Polynomial pSigma;
            const auto xPoly = Polynomial{{0u, 1u}, 1};



/*
 *      r = t/2
 *
 *      sigma0 sigma1 sigma2 sigma3 ...  sigmaR-1 omega0 omega1 ... omegar-1
 *
 *
 */

            Polynomial bigModule;
            bigModule.set(t * 3, 1u);
            auto goppaQuad = goppaPoly.pow(2, bigModule);

#if 1

            Matrix<2 * t, 2 * t, typename GPoly::Member> equationMat;
            Matrix<1, 2 * t, typename GPoly::Member> freeMembers;

            auto syndrome = H2_ * y.transpose();


            for(uint32_t i = 0; i < 2 * t; ++i) {
                typename Polynomial::Member elem;
                for (uint32_t j = 0; j < m; ++j) {
                    elem.set(j, syndrome.at(0, i * m + j));
                }

                pSyndrome.set(i, elem);
            }

            for(uint32_t i = 0; i < t; ++i) {
                Polynomial shift;
                shift.set(i, 1u);
                auto mult = pSyndrome.multiply(shift, goppaQuad);
                for (uint32_t j = 0; j < 2 * t; ++j) {
                    equationMat.at(i, j) = mult.g(j);
                }
            }




            for(uint32_t i = t; i < 2 * t; ++i) {
                for (uint32_t j = 0; j < t; ++j) {
                    if(i == j + t)
                        equationMat.at(i, j) = 1u;
                }
            }
            {
                Polynomial shift;
                shift.set(t, 1u);
                auto mult = pSyndrome.multiply(shift, goppaQuad);
                for (uint32_t i = 0; i < 2 * t; ++i) {
                    freeMembers.at(0, i) = mult.g(i);
                }
            }


            auto almostSigma = inverse(equationMat) * freeMembers;

            for(uint32_t i = 0; i < t; ++i){
                pSigma.set(i, almostSigma.at(0, i));
            }
            pSigma.set(t, 1u);
#else

            auto syndrome = H_ * y.transpose();

            for(uint32_t i = 0; i < t; ++i) {
                typename Polynomial::Member elem;
                for (uint32_t j = 0; j < m; ++j) {
                    elem.set(j, syndrome.at(0, i * m + j));
                }

                pSyndrome.set(i, elem);
            }
            auto pH = pSyndrome.inverse(goppaPoly);

            assert((pH.multiply(pSyndrome, goppaPoly) == Polynomial{{1u}, 1}));

            if(pH == xPoly){
                pSigma = pH;
            } else{
                pH += xPoly;
                auto pD = pH;

                for(uint32_t i = 0; i < m * t - 1u; ++i){
                    pD = pD.pow(2u , goppaPoly);
                }

                assert(pD.pow(2, goppaPoly) == pH && "d^2 == pH must be true");
                auto pInvD = pD.inverse(goppaPoly);
                Polynomial pA, pB;
                pA = pD;
                pB.set(0, 1);
                if(pB.degree() >= pA.degree()) {
                    do {
                        pA.clear();
                        for (uint32_t i = 0; i < t; ++i) {
                            pA.set(i, std::rand() % (pow_const(2u, m) - 1) + 1);

                        }
                        pB = pA.multiply(pInvD, goppaPoly);
                    } while (pB.degree() >= pA.degree());
                }
                assert(pD.multiply(pB, goppaPoly) == pA && "d * b == a must be true");
                pSigma = pA.pow(2u, goppaQuad) + (pB.pow(2u, goppaQuad).multiply(xPoly, goppaQuad));
            }
#endif
            EncryptedMessage errorVec;

            for(uint32_t i = 0; i < n; ++i){
                errorVec.at(i, 0) = pSigma(alphas_.at(i)).isNull() ? 1u : 0u;
            }

            y -= errorVec;
            Matrix<k + 1, n, GaloisElem<2u>> reduceMat;
            reduceMat.copySubmatrix(g_mat_.transpose(), 0, 0);
            reduceMat.template copySubmatrix(y.transpose(), k, 0);


            Message decrypted = reduceMat.reduce().template submatrix<1, k>(k, 0).transpose();

            decrypted = decrypted * inverse(S_);

            return decrypted;
        }

    };



    template<uint32_t n, uint32_t k>
    class McElieceCipher{
        using Key = std::pair<Matrix<n, k, GaloisElem<2u>>, uint32_t>;
        Key pubKey_;
    public:
        using EncryptedMessage = Matrix<n, 1, GaloisElem<2u>>;
        using Message = Matrix<k, 1, GaloisElem<2u>>;
        McElieceCipher(Key pubKey): pubKey_(pubKey){}

        EncryptedMessage encrypt(Message const& m){
            EncryptedMessage z;
            std::deque<uint32_t> freeIds;
            for(uint32_t i = 0; i < n; ++i){
                freeIds.emplace_back(i);
            }


            for(uint32_t i = 0; i < pubKey_.second; ++i){
                uint32_t randPick = randomEngine.get() % freeIds.size();
                z.at(freeIds.at(randPick), 0) = 1u;
                freeIds.erase(freeIds.begin() + randPick);
            }

            return m * pubKey_.first + z;
        }
    };
}

#endif //MCELIECECRYPTO_GALOISFIELD_H
