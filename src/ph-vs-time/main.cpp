#include "y/program.hpp"
#include "y/stream/libc/input.hpp"
#include "y/stream/libc/output.hpp"
#include "y/string.hpp"
#include "y/sequence/vector.hpp"
#include "y/container/algo/crop.hpp"

using namespace Yttrium;

typedef Vector<String,Memory::Pooled> Strings;

class Line : public Object, public Strings
{
public:
    static const char Sep = ';';
    typedef CxxListOf<Line> List;

    explicit Line(const String &input) : Strings(), next(0), prev(0)
    {
        setup(input);
    }

    virtual ~Line() noexcept
    {

    }

    Line *next;
    Line *prev;
private:
    Y_DISABLE_COPY_AND_ASSIGN(Line);
    void setup(const String &input)
    {
        const char * curr = strchr(input.c_str(),Sep);
        if(!curr) return;
        size_t indx=0;
        ++curr;
        while(true)
        {
            const char * const csep = strchr(curr,Sep);
            if(!csep) break;
            String field(curr,csep-curr);
            Algo::Crop(field,isspace);
            for(size_t i=field.size();i>0;--i)
            {
                field[i] = tolower(field[i]);
            }
            pushTail(field);
            curr=csep+1;
            ++indx;
        }
        assert(size()==indx);
        //std::cerr << *this << std::endl;

    }

};


class Lines : public Line::List
{
public:
    inline explicit Lines() noexcept : Line::List() {}
    inline virtual ~Lines() noexcept {}

    inline Lines & operator<<(const String &input)
    {
        pushTail( new Line(input) );
        return *this;
    }

    template <typename FILENAME>
    void loadFrom(const FILENAME &fileName)
    {
        InputFile fp(fileName);
        String    line;
        while(fp.gets(line))
            (*this) << line;

    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(Lines);
};

typedef Vector<double,Memory::Pooled> ColumnType;

class Column : public ColumnType
{
public:
    explicit Column() noexcept : ColumnType() {}
    virtual ~Column() noexcept {}

private:
    Y_DISABLE_COPY_AND_ASSIGN(Column);
};

class Experiment : public Object
{
public:
    typedef CxxListOf<Experiment> List;

    inline explicit Experiment() noexcept : t(), pH(), loading(0), indx(0),  next(0), prev(0) {}
    inline virtual ~Experiment() noexcept {}

    void extract(const Lines &lines, const size_t tidx);
    void generate() const;


    void save() const
    {
        const String fileName = Formatted::Get("xp%u.dat", unsigned(indx));
        OutputFile fp(fileName);
        for(size_t i=1;i<=t.size();++i)
        {
            fp("%.15g %.15g\n", t[i]-loading, pH[i]);
        }
    }


    Column       t;
    Column       pH;
    double       loading;
    const size_t indx;
    Experiment * next;
    Experiment * prev;

private:
    Y_DISABLE_COPY_AND_ASSIGN(Experiment);
};

#include "y/text/ascii/convert.hpp"

void Experiment:: extract(const Lines &lines, const size_t tidx)
{
    assert(tidx>0);
    Coerce(indx)       = tidx;
    const size_t pidx  = tidx+1;
    size_t       iline = 0;
    for(const Line *line=lines.head;line;line=line->next)
    {
        const Readable<String> &arr = *line;
        if(arr.size()<pidx) throw Exception("not enough data for experiment @colum %u", unsigned(tidx));

        const String tString = arr[tidx];
        const String pString = arr[pidx];
        if(tString.size()<=0 || pString.size()<=0 ) continue;
        ++iline;

        std::cerr << "using [" << tString << ";" << pString << "]" << " @" << iline << std::endl;

        if( iline<=1 )
        {
            std::cerr << "skip title" << std::endl;
            continue; // title
        }

        const double tt = ASCII::Convert::ToReal<double>(tString,"time value");

        if( "loading" == pString )
        {
            std::cerr << "found loading@" << tt << std::endl;
            loading = tt;
            continue;
        }

        if( "#num!" == pString )
        {
            continue;
        }

        if( "#div/0!" == pString )
        {
            continue;
        }

        const double pp = ASCII::Convert::ToReal<double>(pString,"pH value");

        if(t.size()>0 && tt<=t.tail()) continue;
        t  << tt;
        pH << pp;
    }

    if( t.size() > 0)
    {
        for(size_t i=t.size();i>0;--i)
        {
            t[i] -= loading;
        }

    }

    loading = 0;

}

#include "y/mkl/interpolation/linear.hpp"


void Experiment:: generate() const
{

    assert(t.size()>0);
    const long t_min = static_cast<long>(ceil(t.head()));
    const long t_max = static_cast<long>(floor(t.tail()));
    if(t_max<=t_min) return;
    const String fn = Formatted::Get("xtr%u.dat",unsigned(indx));
    OutputFile   fp(fn);
    MKL::LinearInterpolation<double> linear;
    for(long i=t_min;i<=t_max;++i)
    {
        const double y = linear(i,t,pH);
        fp("%ld %.15g\n", i, y);
    }


}



Y_Program()
{
    if(argc<=1) {
        std::cerr << "usage: " << program << " filename.csv [time-cols]" << std::endl;
        return 1;
    }

    Lines lines;
    lines.loadFrom(argv[1]);

    Experiment::List xps;
    for(int iarg=2;iarg<argc;++iarg)
    {
        const size_t tidx = ASCII::Convert::To<size_t>(argv[iarg],"tidx");
        std::cerr << "Extracting column with time @" << tidx << std::endl;
        Experiment &xp = *xps.pushTail( new Experiment() );
        xp.extract(lines,tidx);
        xp.save();
        xp.generate();
    }




}
Y_End()
