#include <emscripten.h>

/*extern "C"{
    int convert(const char* inputFileContent);
}*/

/* Compile this program by: 
 * $ g++ -O3 BeEM.cpp -o BeEM
 */

const char* docstring=""
"BeEM input.cif\n"
"    convert PDBx/mmCIF format input file 'input.cif' to Best Effort/Minimal\n"
"    PDB files. Output results to *-pdb-bundle*\n"
"\n"
"option:\n"
"    -p=xxxx          prefix of output file.\n"
"                     default is the PDB ID read from the input\n"
"    -seqres={0,1}    whether to convert SEQRES record\n"
"                     0 - (default) do not convert SEQRES\n"
"                     1 - convert SEQRES\n"
"    -dbref={0,1}     whether to convert dbref record\n"
"                     0 - (default) do not convert DBREF\n"
"                     1 - convert DBREF\n"
"    -gzip={0,1}      whether to perform gzip compression\n"
"                     0 - (default) do not perform compression\n"
"                     1 - perform compression if tar and gzip are available\n"
"    -upper={0,1,2}   whether to convert PDB header text to upper case\n"
"                     0 - do not convert to upper case\n"
"                     1 - (default) only convert header text of single PDB\n"
"                         file to upper case; allow lower case header for\n"
"                         Best Effort/Minimal PDB bundle\n"
"                     2 - convert all PDB text to upper case\n"
"    -maxatom=99999   maximum number of atoms in a file. default is 99999.\n"
"                     no limit on number of atoms if maxatom<=0\n"
" -outfmt={0,1,2,3,4} output format\n"
"                     0 - (default) output a single PDB file if possible;\n"
"                         otherwise, output Best Effort/Minimal PDB bundle\n"
"                     1 - always output Best Effort/Minimal PDB bundle\n"
"                     2 - output one chain per PDB file\n"
"                     3 - always output a single PDB file\n"
"                     4 - output FASTA sequence converted from coordinate\n"
"   -chain=A,B        comma seperated list of chains to output\n"
"                     default is to output all chains\n"
"   -idmap={txt,tsv}  format of chain ID mapping file\n"
"                     txt - (default) space-justified text\n"
"                     tsv - tab-delimited tabular text\n"
"   -ccd5={map,trim}  how to handle expanded chemical component ID >3 characters\n"
"                     map  - (default) map the residue name to reserved set of \n"
"                            chemical component IDs: 01 - 99, DRG, INH, LIG\n"
"                     trim - trim the residue name to keep only the first three\n"
"                            characters\n"
;

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cstdlib>
using namespace std;

/* StringTools START */
string Upper(const string &inputString)
{
    string result=inputString;
    transform(result.begin(), result.end(), result.begin(), ::toupper);
    return result;
}

string Lower(const string &inputString)
{
    string result=inputString;
    transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}

string Join(const string sep, const vector<string>& string_vec,
    const int joinFrom=0)
{
    if (string_vec.size()<=joinFrom) return "";
    string joined_str=string_vec[joinFrom];
    for (int s=joinFrom+1;s<string_vec.size();s++)
        joined_str+=sep+string_vec[s];
    return joined_str;
}

/* split a long string into vectors by whitespace 
 * line          - input string
 * line_vec      - output vector 
 * delimiter     - delimiter */
void Split(const string &line, vector<string> &line_vec,
    const char delimiter=' ',const bool ignore_quotation=false)
{
    bool within_word = false;
    bool within_quotation = false;
    for (size_t pos=0;pos<line.size();pos++)
    {
        if (ignore_quotation==false && (line[pos]=='"' || line[pos]=='\''))
        {
            if (within_quotation) within_quotation=false;
            else within_quotation=true;
        }
        else if (line[pos]=='\n' || line[pos]=='\r')
        {
            within_quotation=false;
        }
        if (line[pos]==delimiter && within_quotation==false)
        {
            within_word = false;
            continue;
        }
        if (!within_word)
        {
            within_word = true;
            line_vec.push_back("");
        }
        line_vec.back()+=line[pos];
    }
}

void clear_line_vec(vector<string> &line_vec)
{
    int i;
    for (i=0;i<line_vec.size();i++) line_vec[i].clear();
    line_vec.clear();
}

string Basename(const string &inputString)
{
    string result="";
    int i;
    vector<string> line_vec;
    Split(inputString,line_vec,'/');
    if (line_vec.size())
    {
        result=line_vec.back();
        clear_line_vec(line_vec);
        Split(result,line_vec,'\\');
        result=line_vec.back();
        clear_line_vec(line_vec);
    }
    return result;
}

string Trim(const string &inputString,const string &char_list=" \n\r\t")
{
    string result = inputString;
    int idxBegin = inputString.find_first_not_of(char_list);
    int idxEnd = inputString.find_last_not_of(char_list);
    if (idxBegin >= 0 && idxEnd >= 0)
        result = inputString.substr(idxBegin, idxEnd + 1 - idxBegin);
    else result = "";
    return result;
}

string lstrip(const string &inputString,const string &char_list=" \n\r\t")
{
    string result = inputString;
    int idxBegin = inputString.find_first_not_of(char_list);
    if (idxBegin >= 0) result = inputString.substr(idxBegin);
    else result = "";
    return result;
}

string rstrip(const string &inputString,const string &char_list=" \n\r\t")
{
    string result=inputString;
    int idxEnd = inputString.find_last_not_of(char_list);
    if (idxEnd >= 0) result = inputString.substr(0, idxEnd + 1);
    else result = "";
    return result;
}

inline bool StartsWith(const string &longString, const string &shortString)
{
    return (longString.size()>=shortString.size() &&
            longString.substr(0,shortString.size())==shortString);
}

inline bool EndsWith(const string &longString, const string &shortString)
{
    return (longString.size()>=shortString.size() &&
            longString.substr(longString.size()-shortString.size(),
                shortString.size())==shortString);
}

inline string formatString(const string &inputString,const int width=8, 
    const int digit=3)
{
    string result=Trim(inputString," ");
    if (StartsWith(result,"00")) result='0'+lstrip(result,"0");
    size_t found=result.find_first_of('.');
    int i;
    if (found==string::npos)
    {
        result+='.';
        found=result.find_first_of('.');
    }
    int curWidth=result.size();
    if (curWidth<found+digit+1)
        for (i=0;i<((found+digit+1)-curWidth);i++) result+='0';
    else if (curWidth>found+digit+1)
    {
        long int extra_prod=1;
        for (i=0;i+found+1<curWidth;i++) extra_prod*=10;
        long int first=atoi((result.substr(0,found)).c_str())*extra_prod;
        long int second=atoi((lstrip(result.substr(found+1),"0")).c_str());
        if (result[0]=='-') second=-second;
        stringstream buf;
        buf<<fixed<<setprecision(digit)<<(first+second+.5)/extra_prod;
        result=buf.str();
        buf.str(string());
    }
    // -0.000
    if (StartsWith(result,"-0."))
    {
        bool allzero=true;
        for (i=found+1;i<result.size();i++)
        {
            if (result[i]!='0')
            {
                allzero=false;
                break;
            }
        }
        if (allzero) result=result.substr(1);
    }
    if (width)
    {
        curWidth=result.size();
        if (curWidth>width)
        {
            result=result.substr(0,width);
            //result=result.substr(result.size()-width);
        }
        else if (curWidth<width)
            for (i=0;i<width-curWidth;i++) result=' '+result;
    }
    return result;
}

/* StringTools END */
/* pstream START */

// PStreams - POSIX Process I/O for C++

//        Copyright (C) 2001 - 2017 Jonathan Wakely
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)
//

/* do not compile on windows, which does not have cygwin */
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) && !defined(__CYGWIN__)
#define NO_PSTREAM
#else

#ifndef REDI_PSTREAM_H_SEEN
#define REDI_PSTREAM_H_SEEN

#include <ios>
#include <streambuf>
#include <istream>
#include <ostream>
#include <string>
#include <vector>
#include <sstream>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>    // for min()
#include <cerrno>       // for errno
#include <cstddef>      // for size_t, NULL
#include <cstdlib>      // for exit()
#include <sys/types.h>  // for pid_t
#include <sys/wait.h>   // for waitpid()
#include <sys/ioctl.h>  // for ioctl() and FIONREAD
#if defined(__sun)
# include <sys/filio.h> // for FIONREAD on Solaris 2.5
#endif
#include <unistd.h>     // for pipe() fork() exec() and filedes functions
#include <signal.h>     // for kill()
#include <fcntl.h>      // for fcntl()
#if REDI_EVISCERATE_PSTREAMS
# include <stdio.h>     // for FILE, fdopen()
#endif


/// The library version.
#define PSTREAMS_VERSION 0x0101   // 1.0.1

/**
 *  @namespace redi
 *  @brief  All PStreams classes are declared in namespace redi.
 *
 *  Like the standard iostreams, PStreams is a set of class templates,
 *  taking a character type and traits type. As with the standard streams
 *  they are most likely to be used with @c char and the default
 *  traits type, so typedefs for this most common case are provided.
 *
 *  The @c pstream_common class template is not intended to be used directly,
 *  it is used internally to provide the common functionality for the
 *  other stream classes.
 */
namespace redi
{
  /// Common base class providing constants and typenames.
  struct pstreams
  {
    /// Type used to specify how to connect to the process.
    typedef std::ios_base::openmode           pmode;

    /// Type used to hold the arguments for a command.
    typedef std::vector<std::string>          argv_type;

    /// Type used for file descriptors.
    typedef int                               fd_type;

    static const pmode pstdin  = std::ios_base::out; ///< Write to stdin
    static const pmode pstdout = std::ios_base::in;  ///< Read from stdout
    static const pmode pstderr = std::ios_base::app; ///< Read from stderr

    /// Create a new process group for the child process.
    static const pmode newpg   = std::ios_base::trunc;

  protected:
    enum { bufsz = 32 };  ///< Size of pstreambuf buffers.
    enum { pbsz  = 2 };   ///< Number of putback characters kept.
  };

  /// Class template for stream buffer.
  template <typename CharT, typename Traits = std::char_traits<CharT> >
    class basic_pstreambuf
    : public std::basic_streambuf<CharT, Traits>
    , public pstreams
    {
    public:
      // Type definitions for dependent types
      typedef CharT                             char_type;
      typedef Traits                            traits_type;
      typedef typename traits_type::int_type    int_type;
      typedef typename traits_type::off_type    off_type;
      typedef typename traits_type::pos_type    pos_type;
      /** @deprecated use pstreams::fd_type instead. */
      typedef fd_type                           fd_t;

      /// Default constructor.
      basic_pstreambuf();

      /// Constructor that initialises the buffer with @a cmd.
      basic_pstreambuf(const std::string& cmd, pmode mode);

      /// Constructor that initialises the buffer with @a file and @a argv.
      basic_pstreambuf( const std::string& file,
                        const argv_type& argv,
                        pmode mode );

      /// Destructor.
      ~basic_pstreambuf();

      /// Initialise the stream buffer with @a cmd.
      basic_pstreambuf*
      open(const std::string& cmd, pmode mode);

      /// Initialise the stream buffer with @a file and @a argv.
      basic_pstreambuf*
      open(const std::string& file, const argv_type& argv, pmode mode);

      /// Close the stream buffer and wait for the process to exit.
      basic_pstreambuf*
      close();

      /// Send a signal to the process.
      basic_pstreambuf*
      kill(int signal = SIGTERM);

      /// Send a signal to the process' process group.
      basic_pstreambuf*
      killpg(int signal = SIGTERM);

      /// Close the pipe connected to the process' stdin.
      void
      peof();

      /// Change active input source.
      bool
      read_err(bool readerr = true);

      /// Report whether the stream buffer has been initialised.
      bool
      is_open() const;

      /// Report whether the process has exited.
      bool
      exited();

#if REDI_EVISCERATE_PSTREAMS
      /// Obtain FILE pointers for each of the process' standard streams.
      std::size_t
      fopen(FILE*& in, FILE*& out, FILE*& err);
#endif

      /// Return the exit status of the process.
      int
      status() const;

      /// Return the error number (errno) for the most recent failed operation.
      int
      error() const;

    protected:
      /// Transfer characters to the pipe when character buffer overflows.
      int_type
      overflow(int_type c);

      /// Transfer characters from the pipe when the character buffer is empty.
      int_type
      underflow();

      /// Make a character available to be returned by the next extraction.
      int_type
      pbackfail(int_type c = traits_type::eof());

      /// Write any buffered characters to the stream.
      int
      sync();

      /// Insert multiple characters into the pipe.
      std::streamsize
      xsputn(const char_type* s, std::streamsize n);

      /// Insert a sequence of characters into the pipe.
      std::streamsize
      write(const char_type* s, std::streamsize n);

      /// Extract a sequence of characters from the pipe.
      std::streamsize
      read(char_type* s, std::streamsize n);

      /// Report how many characters can be read from active input without blocking.
      std::streamsize
      showmanyc();

    protected:
      /// Enumerated type to indicate whether stdout or stderr is to be read.
      enum buf_read_src { rsrc_out = 0, rsrc_err = 1 };

      /// Initialise pipes and fork process.
      pid_t
      fork(pmode mode);

      /// Wait for the child process to exit.
      int
      wait(bool nohang = false);

      /// Return the file descriptor for the output pipe.
      fd_type&
      wpipe();

      /// Return the file descriptor for the active input pipe.
      fd_type&
      rpipe();

      /// Return the file descriptor for the specified input pipe.
      fd_type&
      rpipe(buf_read_src which);

      void
      create_buffers(pmode mode);

      void
      destroy_buffers(pmode mode);

      /// Writes buffered characters to the process' stdin pipe.
      bool
      empty_buffer();

      bool
      fill_buffer(bool non_blocking = false);

      /// Return the active input buffer.
      char_type*
      rbuffer();

      buf_read_src
      switch_read_buffer(buf_read_src);

    private:
      basic_pstreambuf(const basic_pstreambuf&);
      basic_pstreambuf& operator=(const basic_pstreambuf&);

      void
      init_rbuffers();

      pid_t         ppid_;        // pid of process
      fd_type       wpipe_;       // pipe used to write to process' stdin
      fd_type       rpipe_[2];    // two pipes to read from, stdout and stderr
      char_type*    wbuffer_;
      char_type*    rbuffer_[2];
      char_type*    rbufstate_[3];
      /// Index into rpipe_[] to indicate active source for read operations.
      buf_read_src  rsrc_;
      int           status_;      // hold exit status of child process
      int           error_;       // hold errno if fork() or exec() fails
    };

  /// Class template for common base class.
  template <typename CharT, typename Traits = std::char_traits<CharT> >
    class pstream_common
    : virtual public std::basic_ios<CharT, Traits>
    , virtual public pstreams
    {
    protected:
      typedef basic_pstreambuf<CharT, Traits>       streambuf_type;

      typedef pstreams::pmode                       pmode;
      typedef pstreams::argv_type                   argv_type;

      /// Default constructor.
      pstream_common();

      /// Constructor that initialises the stream by starting a process.
      pstream_common(const std::string& cmd, pmode mode);

      /// Constructor that initialises the stream by starting a process.
      pstream_common(const std::string& file, const argv_type& argv, pmode mode);

      /// Pure virtual destructor.
      virtual
      ~pstream_common() = 0;

      /// Start a process.
      void
      do_open(const std::string& cmd, pmode mode);

      /// Start a process.
      void
      do_open(const std::string& file, const argv_type& argv, pmode mode);

    public:
      /// Close the pipe.
      void
      close();

      /// Report whether the stream's buffer has been initialised.
      bool
      is_open() const;

      /// Return the command used to initialise the stream.
      const std::string&
      command() const;

      /// Return a pointer to the stream buffer.
      streambuf_type*
      rdbuf() const;

#if REDI_EVISCERATE_PSTREAMS
      /// Obtain FILE pointers for each of the process' standard streams.
      std::size_t
      fopen(FILE*& in, FILE*& out, FILE*& err);
#endif

    protected:
      std::string       command_; ///< The command used to start the process.
      streambuf_type    buf_;     ///< The stream buffer.
    };


  /**
   * @class basic_ipstream
   * @brief Class template for Input PStreams.
   *
   * Reading from an ipstream reads the command's standard output and/or
   * standard error (depending on how the ipstream is opened)
   * and the command's standard input is the same as that of the process
   * that created the object, unless altered by the command itself.
   */

  template <typename CharT, typename Traits = std::char_traits<CharT> >
    class basic_ipstream
    : public std::basic_istream<CharT, Traits>
    , public pstream_common<CharT, Traits>
    , virtual public pstreams
    {
      typedef std::basic_istream<CharT, Traits>     istream_type;
      typedef pstream_common<CharT, Traits>         pbase_type;

      using pbase_type::buf_;  // declare name in this scope

      // Ensure a basic_ipstream will read from at least one pipe
      pmode readable(pmode mode)
      {
        if (!(mode & (pstdout|pstderr)))
          mode |= pstdout;
        return mode;
      }

    public:
      /// Type used to specify how to connect to the process.
      typedef typename pbase_type::pmode            pmode;

      /// Type used to hold the arguments for a command.
      typedef typename pbase_type::argv_type        argv_type;

      /// Default constructor, creates an uninitialised stream.
      basic_ipstream()
      : istream_type(NULL), pbase_type()
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling do_open() with the supplied
       * arguments.
       *
       * @param cmd   a string containing a shell command.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, pmode)
       */
      explicit
      basic_ipstream(const std::string& cmd, pmode mode = pstdout)
      : istream_type(NULL), pbase_type(cmd, readable(mode))
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling do_open() with the supplied
       * arguments.
       *
       * @param file  a string containing the pathname of a program to execute.
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      basic_ipstream( const std::string& file,
                      const argv_type& argv,
                      pmode mode = pstdout )
      : istream_type(NULL), pbase_type(file, argv, readable(mode))
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling
       * @c do_open(argv[0],argv,mode|pstdout)
       *
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      explicit
      basic_ipstream(const argv_type& argv, pmode mode = pstdout)
      : istream_type(NULL), pbase_type(argv.at(0), argv, readable(mode))
      { }

#if __cplusplus >= 201103L
      template<typename T>
        explicit
        basic_ipstream(std::initializer_list<T> args, pmode mode = pstdout)
        : basic_ipstream(argv_type(args.begin(), args.end()), mode)
        { }
#endif

      /**
       * @brief Destructor.
       *
       * Closes the stream and waits for the child to exit.
       */
      ~basic_ipstream()
      { }

      /**
       * @brief Start a process.
       *
       * Calls do_open( @a cmd , @a mode|pstdout ).
       *
       * @param cmd   a string containing a shell command.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, pmode)
       */
      void
      open(const std::string& cmd, pmode mode = pstdout)
      {
        this->do_open(cmd, readable(mode));
      }

      /**
       * @brief Start a process.
       *
       * Calls do_open( @a file , @a argv , @a mode|pstdout ).
       *
       * @param file  a string containing the pathname of a program to execute.
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      void
      open( const std::string& file,
            const argv_type& argv,
            pmode mode = pstdout )
      {
        this->do_open(file, argv, readable(mode));
      }

      /**
       * @brief Set streambuf to read from process' @c stdout.
       * @return  @c *this
       */
      basic_ipstream&
      out()
      {
        this->buf_.read_err(false);
        return *this;
      }

    };


  /**
   * @class basic_opstream
   * @brief Class template for Output PStreams.
   *
   * Writing to an open opstream writes to the standard input of the command;
   * the command's standard output is the same as that of the process that
   * created the pstream object, unless altered by the command itself.
   */

  template <typename CharT, typename Traits = std::char_traits<CharT> >
    class basic_opstream
    : public std::basic_ostream<CharT, Traits>
    , public pstream_common<CharT, Traits>
    , virtual public pstreams
    {
      typedef std::basic_ostream<CharT, Traits>     ostream_type;
      typedef pstream_common<CharT, Traits>         pbase_type;

      using pbase_type::buf_;  // declare name in this scope

    public:
      /// Type used to specify how to connect to the process.
      typedef typename pbase_type::pmode            pmode;

      /// Type used to hold the arguments for a command.
      typedef typename pbase_type::argv_type        argv_type;

      /// Default constructor, creates an uninitialised stream.
      basic_opstream()
      : ostream_type(NULL), pbase_type()
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling do_open() with the supplied
       * arguments.
       *
       * @param cmd   a string containing a shell command.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, pmode)
       */
      explicit
      basic_opstream(const std::string& cmd, pmode mode = pstdin)
      : ostream_type(NULL), pbase_type(cmd, mode|pstdin)
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling do_open() with the supplied
       * arguments.
       *
       * @param file  a string containing the pathname of a program to execute.
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      basic_opstream( const std::string& file,
                      const argv_type& argv,
                      pmode mode = pstdin )
      : ostream_type(NULL), pbase_type(file, argv, mode|pstdin)
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling
       * @c do_open(argv[0],argv,mode|pstdin)
       *
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      explicit
      basic_opstream(const argv_type& argv, pmode mode = pstdin)
      : ostream_type(NULL), pbase_type(argv.at(0), argv, mode|pstdin)
      { }

#if __cplusplus >= 201103L
      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * @param args  a list of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      template<typename T>
        explicit
        basic_opstream(std::initializer_list<T> args, pmode mode = pstdin)
        : basic_opstream(argv_type(args.begin(), args.end()), mode)
        { }
#endif

      /**
       * @brief Destructor
       *
       * Closes the stream and waits for the child to exit.
       */
      ~basic_opstream() { }

      /**
       * @brief Start a process.
       *
       * Calls do_open( @a cmd , @a mode|pstdin ).
       *
       * @param cmd   a string containing a shell command.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, pmode)
       */
      void
      open(const std::string& cmd, pmode mode = pstdin)
      {
        this->do_open(cmd, mode|pstdin);
      }

      /**
       * @brief Start a process.
       *
       * Calls do_open( @a file , @a argv , @a mode|pstdin ).
       *
       * @param file  a string containing the pathname of a program to execute.
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      void
      open( const std::string& file,
            const argv_type& argv,
            pmode mode = pstdin)
      {
        this->do_open(file, argv, mode|pstdin);
      }
    };


  /**
   * @class basic_pstream
   * @brief Class template for Bidirectional PStreams.
   *
   * Writing to a pstream opened with @c pmode @c pstdin writes to the
   * standard input of the command.
   * Reading from a pstream opened with @c pmode @c pstdout and/or @c pstderr
   * reads the command's standard output and/or standard error.
   * Any of the process' @c stdin, @c stdout or @c stderr that is not
   * connected to the pstream (as specified by the @c pmode)
   * will be the same as the process that created the pstream object,
   * unless altered by the command itself.
   */
  template <typename CharT, typename Traits = std::char_traits<CharT> >
    class basic_pstream
    : public std::basic_iostream<CharT, Traits>
    , public pstream_common<CharT, Traits>
    , virtual public pstreams
    {
      typedef std::basic_iostream<CharT, Traits>    iostream_type;
      typedef pstream_common<CharT, Traits>         pbase_type;

      using pbase_type::buf_;  // declare name in this scope

    public:
      /// Type used to specify how to connect to the process.
      typedef typename pbase_type::pmode            pmode;

      /// Type used to hold the arguments for a command.
      typedef typename pbase_type::argv_type        argv_type;

      /// Default constructor, creates an uninitialised stream.
      basic_pstream()
      : iostream_type(NULL), pbase_type()
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling do_open() with the supplied
       * arguments.
       *
       * @param cmd   a string containing a shell command.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, pmode)
       */
      explicit
      basic_pstream(const std::string& cmd, pmode mode = pstdout|pstdin)
      : iostream_type(NULL), pbase_type(cmd, mode)
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling do_open() with the supplied
       * arguments.
       *
       * @param file  a string containing the pathname of a program to execute.
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      basic_pstream( const std::string& file,
                     const argv_type& argv,
                     pmode mode = pstdout|pstdin )
      : iostream_type(NULL), pbase_type(file, argv, mode)
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling
       * @c do_open(argv[0],argv,mode)
       *
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      explicit
      basic_pstream(const argv_type& argv, pmode mode = pstdout|pstdin)
      : iostream_type(NULL), pbase_type(argv.at(0), argv, mode)
      { }

#if __cplusplus >= 201103L
      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * @param l     a list of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      template<typename T>
        explicit
        basic_pstream(std::initializer_list<T> l, pmode mode = pstdout|pstdin)
        : basic_pstream(argv_type(l.begin(), l.end()), mode)
        { }
#endif

      /**
       * @brief Destructor
       *
       * Closes the stream and waits for the child to exit.
       */
      ~basic_pstream() { }

      /**
       * @brief Start a process.
       *
       * Calls do_open( @a cnd , @a mode ).
       *
       * @param cmd   a string containing a shell command.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, pmode)
       */
      void
      open(const std::string& cmd, pmode mode = pstdout|pstdin)
      {
        this->do_open(cmd, mode);
      }

      /**
       * @brief Start a process.
       *
       * Calls do_open( @a file , @a argv , @a mode ).
       *
       * @param file  a string containing the pathname of a program to execute.
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      void
      open( const std::string& file,
            const argv_type& argv,
            pmode mode = pstdout|pstdin )
      {
        this->do_open(file, argv, mode);
      }

      /**
       * @brief Set streambuf to read from process' @c stdout.
       * @return  @c *this
       */
      basic_pstream&
      out()
      {
        this->buf_.read_err(false);
        return *this;
      }
    };


  /**
   * @class basic_rpstream
   * @brief Class template for Restricted PStreams.
   *
   * Writing to an rpstream opened with @c pmode @c pstdin writes to the
   * standard input of the command.
   * It is not possible to read directly from an rpstream object, to use
   * an rpstream as in istream you must call either basic_rpstream::out()
   * or basic_rpstream::err(). This is to prevent accidental reads from
   * the wrong input source. If the rpstream was not opened with @c pmode
   * @c pstderr then the class cannot read the process' @c stderr, and
   * basic_rpstream::err() will return an istream that reads from the
   * process' @c stdout, and vice versa.
   * Reading from an rpstream opened with @c pmode @c pstdout and/or
   * @c pstderr reads the command's standard output and/or standard error.
   * Any of the process' @c stdin, @c stdout or @c stderr that is not
   * connected to the pstream (as specified by the @c pmode)
   * will be the same as the process that created the pstream object,
   * unless altered by the command itself.
   */

  template <typename CharT, typename Traits = std::char_traits<CharT> >
    class basic_rpstream
    : public std::basic_ostream<CharT, Traits>
    , private std::basic_istream<CharT, Traits>
    , private pstream_common<CharT, Traits>
    , virtual public pstreams
    {
      typedef std::basic_ostream<CharT, Traits>     ostream_type;
      typedef std::basic_istream<CharT, Traits>     istream_type;
      typedef pstream_common<CharT, Traits>         pbase_type;

      using pbase_type::buf_;  // declare name in this scope

    public:
      /// Type used to specify how to connect to the process.
      typedef typename pbase_type::pmode            pmode;

      /// Type used to hold the arguments for a command.
      typedef typename pbase_type::argv_type        argv_type;

      /// Default constructor, creates an uninitialised stream.
      basic_rpstream()
      : ostream_type(NULL), istream_type(NULL), pbase_type()
      { }

      /**
       * @brief  Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling do_open() with the supplied
       * arguments.
       *
       * @param cmd   a string containing a shell command.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, pmode)
       */
      explicit
      basic_rpstream(const std::string& cmd, pmode mode = pstdout|pstdin)
      : ostream_type(NULL) , istream_type(NULL) , pbase_type(cmd, mode)
      { }

      /**
       * @brief  Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling do_open() with the supplied
       * arguments.
       *
       * @param file a string containing the pathname of a program to execute.
       * @param argv a vector of argument strings passed to the new program.
       * @param mode the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      basic_rpstream( const std::string& file,
                      const argv_type& argv,
                      pmode mode = pstdout|pstdin )
      : ostream_type(NULL), istream_type(NULL), pbase_type(file, argv, mode)
      { }

      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * Initialises the stream buffer by calling
       * @c do_open(argv[0],argv,mode)
       *
       * @param argv  a vector of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      explicit
      basic_rpstream(const argv_type& argv, pmode mode = pstdout|pstdin)
      : ostream_type(NULL), istream_type(NULL),
        pbase_type(argv.at(0), argv, mode)
      { }

#if __cplusplus >= 201103L
      /**
       * @brief Constructor that initialises the stream by starting a process.
       *
       * @param l     a list of argument strings passed to the new program.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      template<typename T>
        explicit
        basic_rpstream(std::initializer_list<T> l, pmode mode = pstdout|pstdin)
        : basic_rpstream(argv_type(l.begin(), l.end()), mode)
        { }
#endif

      /// Destructor
      ~basic_rpstream() { }

      /**
       * @brief  Start a process.
       *
       * Calls do_open( @a cmd , @a mode ).
       *
       * @param cmd   a string containing a shell command.
       * @param mode  the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, pmode)
       */
      void
      open(const std::string& cmd, pmode mode = pstdout|pstdin)
      {
        this->do_open(cmd, mode);
      }

      /**
       * @brief  Start a process.
       *
       * Calls do_open( @a file , @a argv , @a mode ).
       *
       * @param file a string containing the pathname of a program to execute.
       * @param argv a vector of argument strings passed to the new program.
       * @param mode the I/O mode to use when opening the pipe.
       * @see   do_open(const std::string&, const argv_type&, pmode)
       */
      void
      open( const std::string& file,
            const argv_type& argv,
            pmode mode = pstdout|pstdin )
      {
        this->do_open(file, argv, mode);
      }

      /**
       * @brief  Obtain a reference to the istream that reads
       *         the process' @c stdout.
       * @return @c *this
       */
      istream_type&
      out()
      {
        this->buf_.read_err(false);
        return *this;
      }
    };


  /// Type definition for common template specialisation.
  typedef basic_pstreambuf<char> pstreambuf;
  /// Type definition for common template specialisation.
  typedef basic_ipstream<char> ipstream;
  /// Type definition for common template specialisation.
  typedef basic_opstream<char> opstream;
  /// Type definition for common template specialisation.
  typedef basic_pstream<char> pstream;
  /// Type definition for common template specialisation.
  typedef basic_rpstream<char> rpstream;


  /**
   * When inserted into an output pstream the manipulator calls
   * basic_pstreambuf<C,T>::peof() to close the output pipe,
   * causing the child process to receive the end-of-file indicator
   * on subsequent reads from its @c stdin stream.
   *
   * @brief   Manipulator to close the pipe connected to the process' stdin.
   * @param   s  An output PStream class.
   * @return  The stream object the manipulator was invoked on.
   * @see     basic_pstreambuf<C,T>::peof()
   * @relates basic_opstream basic_pstream basic_rpstream
   */
  template <typename C, typename T>
    inline std::basic_ostream<C,T>&
    peof(std::basic_ostream<C,T>& s)
    {
      typedef basic_pstreambuf<C,T> pstreambuf_type;
      if (pstreambuf_type* p = dynamic_cast<pstreambuf_type*>(s.rdbuf()))
        p->peof();
      return s;
    }


  /*
   * member definitions for pstreambuf
   */


  /**
   * @class basic_pstreambuf
   * Provides underlying streambuf functionality for the PStreams classes.
   */

  /** Creates an uninitialised stream buffer. */
  template <typename C, typename T>
    inline
    basic_pstreambuf<C,T>::basic_pstreambuf()
    : ppid_(-1)   // initialise to -1 to indicate no process run yet.
    , wpipe_(-1)
    , wbuffer_(NULL)
    , rsrc_(rsrc_out)
    , status_(-1)
    , error_(0)
    {
      init_rbuffers();
    }

  /**
   * Initialises the stream buffer by calling open() with the supplied
   * arguments.
   *
   * @param cmd   a string containing a shell command.
   * @param mode  the I/O mode to use when opening the pipe.
   * @see   open()
   */
  template <typename C, typename T>
    inline
    basic_pstreambuf<C,T>::basic_pstreambuf(const std::string& cmd, pmode mode)
    : ppid_(-1)   // initialise to -1 to indicate no process run yet.
    , wpipe_(-1)
    , wbuffer_(NULL)
    , rsrc_(rsrc_out)
    , status_(-1)
    , error_(0)
    {
      init_rbuffers();
      open(cmd, mode);
    }

  /**
   * Initialises the stream buffer by calling open() with the supplied
   * arguments.
   *
   * @param file  a string containing the name of a program to execute.
   * @param argv  a vector of argument strings passsed to the new program.
   * @param mode  the I/O mode to use when opening the pipe.
   * @see   open()
   */
  template <typename C, typename T>
    inline
    basic_pstreambuf<C,T>::basic_pstreambuf( const std::string& file,
                                             const argv_type& argv,
                                             pmode mode )
    : ppid_(-1)   // initialise to -1 to indicate no process run yet.
    , wpipe_(-1)
    , wbuffer_(NULL)
    , rsrc_(rsrc_out)
    , status_(-1)
    , error_(0)
    {
      init_rbuffers();
      open(file, argv, mode);
    }

  /**
   * Closes the stream by calling close().
   * @see close()
   */
  template <typename C, typename T>
    inline
    basic_pstreambuf<C,T>::~basic_pstreambuf()
    {
      close();
    }

  /**
   * Starts a new process by passing @a command to the shell (/bin/sh)
   * and opens pipes to the process with the specified @a mode.
   *
   * If @a mode contains @c pstdout the initial read source will be
   * the child process' stdout, otherwise if @a mode  contains @c pstderr
   * the initial read source will be the child's stderr.
   *
   * Will duplicate the actions of  the  shell  in searching for an
   * executable file if the specified file name does not contain a slash (/)
   * character.
   *
   * @warning
   * There is no way to tell whether the shell command succeeded, this
   * function will always succeed unless resource limits (such as
   * memory usage, or number of processes or open files) are exceeded.
   * This means is_open() will return true even if @a command cannot
   * be executed.
   * Use pstreambuf::open(const std::string&, const argv_type&, pmode)
   * if you need to know whether the command failed to execute.
   *
   * @param   command  a string containing a shell command.
   * @param   mode     a bitwise OR of one or more of @c out, @c in, @c err.
   * @return  NULL if the shell could not be started or the
   *          pipes could not be opened, @c this otherwise.
   * @see     <b>execl</b>(3)
   */
  template <typename C, typename T>
    basic_pstreambuf<C,T>*
    basic_pstreambuf<C,T>::open(const std::string& command, pmode mode)
    {
      const char * shell_path = "/bin/sh";
#if 0
      const std::string argv[] = { "sh", "-c", command };
      return this->open(shell_path, argv_type(argv, argv+3), mode);
#else
      basic_pstreambuf<C,T>* ret = NULL;

      if (!is_open())
      {
        switch(fork(mode))
        {
        case 0 :
          // this is the new process, exec command
          ::execl(shell_path, "sh", "-c", command.c_str(), (char*)NULL);

          // can only reach this point if exec() failed

          // parent can get exit code from waitpid()
          ::_exit(errno);
          // using std::exit() would make static dtors run twice

        case -1 :
          // couldn't fork, error already handled in pstreambuf::fork()
          break;

        default :
          // this is the parent process
          // activate buffers
          create_buffers(mode);
          ret = this;
        }
      }
      return ret;
#endif
    }

  /**
   * @brief  Helper function to close a file descriptor.
   *
   * Inspects @a fd and calls <b>close</b>(3) if it has a non-negative value.
   *
   * @param   fd  a file descriptor.
   * @relates basic_pstreambuf
   */
  inline void
  close_fd(pstreams::fd_type& fd)
  {
    if (fd >= 0 && ::close(fd) == 0)
      fd = -1;
  }

  /**
   * @brief  Helper function to close an array of file descriptors.
   *
   * Calls @c close_fd() on each member of the array.
   * The length of the array is determined automatically by
   * template argument deduction to avoid errors.
   *
   * @param   fds  an array of file descriptors.
   * @relates basic_pstreambuf
   */
  template <int N>
    inline void
    close_fd_array(pstreams::fd_type (&fds)[N])
    {
      for (std::size_t i = 0; i < N; ++i)
        close_fd(fds[i]);
    }

  /**
   * Starts a new process by executing @a file with the arguments in
   * @a argv and opens pipes to the process with the specified @a mode.
   *
   * By convention @c argv[0] should be the file name of the file being
   * executed.
   *
   * If @a mode contains @c pstdout the initial read source will be
   * the child process' stdout, otherwise if @a mode  contains @c pstderr
   * the initial read source will be the child's stderr.
   *
   * Will duplicate the actions of  the  shell  in searching for an
   * executable file if the specified file name does not contain a slash (/)
   * character.
   *
   * Iff @a file is successfully executed then is_open() will return true.
   * Otherwise, pstreambuf::error() can be used to obtain the value of
   * @c errno that was set by <b>execvp</b>(3) in the child process.
   *
   * The exit status of the new process will be returned by
   * pstreambuf::status() after pstreambuf::exited() returns true.
   *
   * @param   file  a string containing the pathname of a program to execute.
   * @param   argv  a vector of argument strings passed to the new program.
   * @param   mode  a bitwise OR of one or more of @c out, @c in and @c err.
   * @return  NULL if a pipe could not be opened or if the program could
   *          not be executed, @c this otherwise.
   * @see     <b>execvp</b>(3)
   */
  template <typename C, typename T>
    basic_pstreambuf<C,T>*
    basic_pstreambuf<C,T>::open( const std::string& file,
                                 const argv_type& argv,
                                 pmode mode )
    {
      basic_pstreambuf<C,T>* ret = NULL;

      if (!is_open())
      {
        // constants for read/write ends of pipe
        enum { RD, WR };

        // open another pipe and set close-on-exec
        fd_type ck_exec[] = { -1, -1 };
        if (-1 == ::pipe(ck_exec)
            || -1 == ::fcntl(ck_exec[RD], F_SETFD, FD_CLOEXEC)
            || -1 == ::fcntl(ck_exec[WR], F_SETFD, FD_CLOEXEC))
        {
          error_ = errno;
          close_fd_array(ck_exec);
        }
        else
        {
          switch(fork(mode))
          {
          case 0 :
            // this is the new process, exec command
            {
              char** arg_v = new char*[argv.size()+1];
              for (std::size_t i = 0; i < argv.size(); ++i)
              {
                const std::string& src = argv[i];
                char*& dest = arg_v[i];
                dest = new char[src.size()+1];
                dest[ src.copy(dest, src.size()) ] = '\0';
              }
              arg_v[argv.size()] = NULL;

              ::execvp(file.c_str(), arg_v);

              // can only reach this point if exec() failed

              // parent can get error code from ck_exec pipe
              error_ = errno;

              while (::write(ck_exec[WR], &error_, sizeof(error_)) == -1
                  && errno == EINTR)
              { }

              ::close(ck_exec[WR]);
              ::close(ck_exec[RD]);

              ::_exit(error_);
              // using std::exit() would make static dtors run twice
            }

          case -1 :
            // couldn't fork, error already handled in pstreambuf::fork()
            close_fd_array(ck_exec);
            break;

          default :
            // this is the parent process

            // check child called exec() successfully
            ::close(ck_exec[WR]);
            switch (::read(ck_exec[RD], &error_, sizeof(error_)))
            {
            case 0:
              // activate buffers
              create_buffers(mode);
              ret = this;
              break;
            case -1:
              error_ = errno;
              break;
            default:
              // error_ contains error code from child
              // call wait() to clean up and set ppid_ to 0
              this->wait();
              break;
            }
            ::close(ck_exec[RD]);
          }
        }
      }
      return ret;
    }

  /**
   * Creates pipes as specified by @a mode and calls @c fork() to create
   * a new process. If the fork is successful the parent process stores
   * the child's PID and the opened pipes and the child process replaces
   * its standard streams with the opened pipes.
   *
   * If an error occurs the error code will be set to one of the possible
   * errors for @c pipe() or @c fork().
   * See your system's documentation for these error codes.
   *
   * @param   mode  an OR of pmodes specifying which of the child's
   *                standard streams to connect to.
   * @return  On success the PID of the child is returned in the parent's
   *          context and zero is returned in the child's context.
   *          On error -1 is returned and the error code is set appropriately.
   */
  template <typename C, typename T>
    pid_t
    basic_pstreambuf<C,T>::fork(pmode mode)
    {
      pid_t pid = -1;

      // Three pairs of file descriptors, for pipes connected to the
      // process' stdin, stdout and stderr
      // (stored in a single array so close_fd_array() can close all at once)
      fd_type fd[] = { -1, -1, -1, -1, -1, -1 };
      fd_type* const pin = fd;
      fd_type* const pout = fd+2;
      fd_type* const perr = fd+4;

      // constants for read/write ends of pipe
      enum { RD, WR };

      // N.B.
      // For the pstreambuf pin is an output stream and
      // pout and perr are input streams.

      if (!error_ && mode&pstdin && ::pipe(pin))
        error_ = errno;

      if (!error_ && mode&pstdout && ::pipe(pout))
        error_ = errno;

      if (!error_ && mode&pstderr && ::pipe(perr))
        error_ = errno;

      if (!error_)
      {
        pid = ::fork();
        switch (pid)
        {
          case 0 :
          {
            // this is the new process

            // for each open pipe close one end and redirect the
            // respective standard stream to the other end

            if (*pin >= 0)
            {
              ::close(pin[WR]);
              ::dup2(pin[RD], STDIN_FILENO);
              ::close(pin[RD]);
            }
            if (*pout >= 0)
            {
              ::close(pout[RD]);
              ::dup2(pout[WR], STDOUT_FILENO);
              ::close(pout[WR]);
            }
            if (*perr >= 0)
            {
              ::close(perr[RD]);
              ::dup2(perr[WR], STDERR_FILENO);
              ::close(perr[WR]);
            }

#ifdef _POSIX_JOB_CONTROL
            if (mode&newpg)
              ::setpgid(0, 0); // Change to a new process group
#endif

            break;
          }
          case -1 :
          {
            // couldn't fork for some reason
            error_ = errno;
            // close any open pipes
            close_fd_array(fd);
            break;
          }
          default :
          {
            // this is the parent process, store process' pid
            ppid_ = pid;

            // store one end of open pipes and close other end
            if (*pin >= 0)
            {
              wpipe_ = pin[WR];
              ::close(pin[RD]);
            }
            if (*pout >= 0)
            {
              rpipe_[rsrc_out] = pout[RD];
              ::close(pout[WR]);
            }
            if (*perr >= 0)
            {
              rpipe_[rsrc_err] = perr[RD];
              ::close(perr[WR]);
            }
          }
        }
      }
      else
      {
        // close any pipes we opened before failure
        close_fd_array(fd);
      }
      return pid;
    }

  /**
   * Closes all pipes and calls wait() to wait for the process to finish.
   * If an error occurs the error code will be set to one of the possible
   * errors for @c waitpid().
   * See your system's documentation for these errors.
   *
   * @return  @c this on successful close or @c NULL if there is no
   *          process to close or if an error occurs.
   */
  template <typename C, typename T>
    basic_pstreambuf<C,T>*
    basic_pstreambuf<C,T>::close()
    {
      const bool running = is_open();

      sync(); // this might call wait() and reap the child process

      // rather than trying to work out whether or not we need to clean up
      // just do it anyway, all cleanup functions are safe to call twice.

      destroy_buffers(pstdin|pstdout|pstderr);

      // close pipes before wait() so child gets EOF/SIGPIPE
      close_fd(wpipe_);
      close_fd_array(rpipe_);

      do
      {
        error_ = 0;
      } while (wait() == -1 && error() == EINTR);

      return running ? this : NULL;
    }

  /**
   *  Called on construction to initialise the arrays used for reading.
   */
  template <typename C, typename T>
    inline void
    basic_pstreambuf<C,T>::init_rbuffers()
    {
      rpipe_[rsrc_out] = rpipe_[rsrc_err] = -1;
      rbuffer_[rsrc_out] = rbuffer_[rsrc_err] = NULL;
      rbufstate_[0] = rbufstate_[1] = rbufstate_[2] = NULL;
    }

  template <typename C, typename T>
    void
    basic_pstreambuf<C,T>::create_buffers(pmode mode)
    {
      if (mode & pstdin)
      {
        delete[] wbuffer_;
        wbuffer_ = new char_type[bufsz];
        this->setp(wbuffer_, wbuffer_ + bufsz);
      }
      if (mode & pstdout)
      {
        delete[] rbuffer_[rsrc_out];
        rbuffer_[rsrc_out] = new char_type[bufsz];
        rsrc_ = rsrc_out;
        this->setg(rbuffer_[rsrc_out] + pbsz, rbuffer_[rsrc_out] + pbsz,
            rbuffer_[rsrc_out] + pbsz);
      }
      if (mode & pstderr)
      {
        delete[] rbuffer_[rsrc_err];
        rbuffer_[rsrc_err] = new char_type[bufsz];
        if (!(mode & pstdout))
        {
          rsrc_ = rsrc_err;
          this->setg(rbuffer_[rsrc_err] + pbsz, rbuffer_[rsrc_err] + pbsz,
              rbuffer_[rsrc_err] + pbsz);
        }
      }
    }

  template <typename C, typename T>
    void
    basic_pstreambuf<C,T>::destroy_buffers(pmode mode)
    {
      if (mode & pstdin)
      {
        this->setp(NULL, NULL);
        delete[] wbuffer_;
        wbuffer_ = NULL;
      }
      if (mode & pstdout)
      {
        if (rsrc_ == rsrc_out)
          this->setg(NULL, NULL, NULL);
        delete[] rbuffer_[rsrc_out];
        rbuffer_[rsrc_out] = NULL;
      }
      if (mode & pstderr)
      {
        if (rsrc_ == rsrc_err)
          this->setg(NULL, NULL, NULL);
        delete[] rbuffer_[rsrc_err];
        rbuffer_[rsrc_err] = NULL;
      }
    }

  template <typename C, typename T>
    typename basic_pstreambuf<C,T>::buf_read_src
    basic_pstreambuf<C,T>::switch_read_buffer(buf_read_src src)
    {
      if (rsrc_ != src)
      {
        char_type* tmpbufstate[] = {this->eback(), this->gptr(), this->egptr()};
        this->setg(rbufstate_[0], rbufstate_[1], rbufstate_[2]);
        for (std::size_t i = 0; i < 3; ++i)
          rbufstate_[i] = tmpbufstate[i];
        rsrc_ = src;
      }
      return rsrc_;
    }

  /**
   * Suspends execution and waits for the associated process to exit, or
   * until a signal is delivered whose action is to terminate the current
   * process or to call a signal handling function. If the process has
   * already exited (i.e. it is a "zombie" process) then wait() returns
   * immediately.  Waiting for the child process causes all its system
   * resources to be freed.
   *
   * error() will return EINTR if wait() is interrupted by a signal.
   *
   * @param   nohang  true to return immediately if the process has not exited.
   * @return  1 if the process has exited and wait() has not yet been called.
   *          0 if @a nohang is true and the process has not exited yet.
   *          -1 if no process has been started or if an error occurs,
   *          in which case the error can be found using error().
   */
  template <typename C, typename T>
    int
    basic_pstreambuf<C,T>::wait(bool nohang)
    {
      int child_exited = -1;
      if (is_open())
      {
        int exit_status;
        switch(::waitpid(ppid_, &exit_status, nohang ? WNOHANG : 0))
        {
          case 0 :
            // nohang was true and process has not exited
            child_exited = 0;
            break;
          case -1 :
            error_ = errno;
            break;
          default :
            // process has exited
            ppid_ = 0;
            status_ = exit_status;
            child_exited = 1;
            // Close wpipe, would get SIGPIPE if we used it.
            destroy_buffers(pstdin);
            close_fd(wpipe_);
            // Must free read buffers and pipes on destruction
            // or next call to open()/close()
            break;
        }
      }
      return child_exited;
    }

  /**
   * Sends the specified signal to the process.  A signal can be used to
   * terminate a child process that would not exit otherwise.
   *
   * If an error occurs the error code will be set to one of the possible
   * errors for @c kill().  See your system's documentation for these errors.
   *
   * @param   signal  A signal to send to the child process.
   * @return  @c this or @c NULL if @c kill() fails.
   */
  template <typename C, typename T>
    inline basic_pstreambuf<C,T>*
    basic_pstreambuf<C,T>::kill(int signal)
    {
      basic_pstreambuf<C,T>* ret = NULL;
      if (is_open())
      {
        if (::kill(ppid_, signal))
          error_ = errno;
        else
        {
#if 0
          // TODO call exited() to check for exit and clean up? leave to user?
          if (signal==SIGTERM || signal==SIGKILL)
            this->exited();
#endif
          ret = this;
        }
      }
      return ret;
    }

  /**
   * Sends the specified signal to the process group of the child process.
   * A signal can be used to terminate a child process that would not exit
   * otherwise, or to kill the process and its own children.
   *
   * If an error occurs the error code will be set to one of the possible
   * errors for @c getpgid() or @c kill().  See your system's documentation
   * for these errors. If the child is in the current process group then
   * NULL will be returned and the error code set to EPERM.
   *
   * @param   signal  A signal to send to the child process.
   * @return  @c this on success or @c NULL on failure.
   */
  template <typename C, typename T>
    inline basic_pstreambuf<C,T>*
    basic_pstreambuf<C,T>::killpg(int signal)
    {
      basic_pstreambuf<C,T>* ret = NULL;
#ifdef _POSIX_JOB_CONTROL
      if (is_open())
      {
        pid_t pgid = ::getpgid(ppid_);
        if (pgid == -1)
          error_ = errno;
        else if (pgid == ::getpgrp())
          error_ = EPERM;  // Don't commit suicide
        else if (::killpg(pgid, signal))
          error_ = errno;
        else
          ret = this;
      }
#else
      error_ = ENOTSUP;
#endif
      return ret;
    }

  /**
   *  This function can call pstreambuf::wait() and so may change the
   *  object's state if the child process has already exited.
   *
   *  @return  True if the associated process has exited, false otherwise.
   *  @see     basic_pstreambuf<C,T>::wait()
   */
  template <typename C, typename T>
    inline bool
    basic_pstreambuf<C,T>::exited()
    {
      return ppid_ == 0 || wait(true)==1;
    }

  /**
   *  @return  The error code of the most recently failed operation, or zero.
   */
  template <typename C, typename T>
    inline int
    basic_pstreambuf<C,T>::error() const
    {
      return error_;
    }

  /**
   *  Closes the output pipe, causing the child process to receive the
   *  end-of-file indicator on subsequent reads from its @c stdin stream.
   */
  template <typename C, typename T>
    inline void
    basic_pstreambuf<C,T>::peof()
    {
      sync();
      destroy_buffers(pstdin);
      close_fd(wpipe_);
    }

  /**
   * Unlike pstreambuf::exited(), this function will not call wait() and
   * so will not change the object's state.  This means that once a child
   * process is executed successfully this function will continue to
   * return true even after the process exits (until wait() is called.)
   *
   * @return  true if a previous call to open() succeeded and wait() has
   *          not been called and determined that the process has exited,
   *          false otherwise.
   */
  template <typename C, typename T>
    inline bool
    basic_pstreambuf<C,T>::is_open() const
    {
      return ppid_ > 0;
    }

  /**
   * Toggle the stream used for reading. If @a readerr is @c true then the
   * process' @c stderr output will be used for subsequent extractions, if
   * @a readerr is false the the process' stdout will be used.
   * @param   readerr  @c true to read @c stderr, @c false to read @c stdout.
   * @return  @c true if the requested stream is open and will be used for
   *          subsequent extractions, @c false otherwise.
   */
  template <typename C, typename T>
    inline bool
    basic_pstreambuf<C,T>::read_err(bool readerr)
    {
      buf_read_src src = readerr ? rsrc_err : rsrc_out;
      if (rpipe_[src]>=0)
      {
        switch_read_buffer(src);
        return true;
      }
      return false;
    }

  /**
   * Called when the internal character buffer is not present or is full,
   * to transfer the buffer contents to the pipe.
   *
   * @param   c  a character to be written to the pipe.
   * @return  @c traits_type::eof() if an error occurs, otherwise if @a c
   *          is not equal to @c traits_type::eof() it will be buffered and
   *          a value other than @c traits_type::eof() returned to indicate
   *          success.
   */
  template <typename C, typename T>
    typename basic_pstreambuf<C,T>::int_type
    basic_pstreambuf<C,T>::overflow(int_type c)
    {
      if (!empty_buffer())
        return traits_type::eof();
      else if (!traits_type::eq_int_type(c, traits_type::eof()))
        return this->sputc(c);
      else
        return traits_type::not_eof(c);
    }


  template <typename C, typename T>
    int
    basic_pstreambuf<C,T>::sync()
    {
      return !exited() && empty_buffer() ? 0 : -1;
    }

  /**
   * @param   s  character buffer.
   * @param   n  buffer length.
   * @return  the number of characters written.
   */
  template <typename C, typename T>
    std::streamsize
    basic_pstreambuf<C,T>::xsputn(const char_type* s, std::streamsize n)
    {
      std::streamsize done = 0;
      while (done < n)
      {
        if (std::streamsize nbuf = this->epptr() - this->pptr())
        {
          nbuf = std::min(nbuf, n - done);
          traits_type::copy(this->pptr(), s + done, nbuf);
          this->pbump(nbuf);
          done += nbuf;
        }
        else if (!empty_buffer())
          break;
      }
      return done;
    }

  /**
   * @return  true if the buffer was emptied, false otherwise.
   */
  template <typename C, typename T>
    bool
    basic_pstreambuf<C,T>::empty_buffer()
    {
      const std::streamsize count = this->pptr() - this->pbase();
      if (count > 0)
      {
        const std::streamsize written = this->write(this->wbuffer_, count);
        if (written > 0)
        {
          if (const std::streamsize unwritten = count - written)
            traits_type::move(this->pbase(), this->pbase()+written, unwritten);
          this->pbump(-written);
          return true;
        }
      }
      return false;
    }

  /**
   * Called when the internal character buffer is is empty, to re-fill it
   * from the pipe.
   *
   * @return The first available character in the buffer,
   * or @c traits_type::eof() in case of failure.
   */
  template <typename C, typename T>
    typename basic_pstreambuf<C,T>::int_type
    basic_pstreambuf<C,T>::underflow()
    {
      if (this->gptr() < this->egptr() || fill_buffer())
        return traits_type::to_int_type(*this->gptr());
      else
        return traits_type::eof();
    }

  /**
   * Attempts to make @a c available as the next character to be read by
   * @c sgetc().
   *
   * @param   c   a character to make available for extraction.
   * @return  @a c if the character can be made available,
   *          @c traits_type::eof() otherwise.
   */
  template <typename C, typename T>
    typename basic_pstreambuf<C,T>::int_type
    basic_pstreambuf<C,T>::pbackfail(int_type c)
    {
      if (this->gptr() != this->eback())
      {
        this->gbump(-1);
        if (!traits_type::eq_int_type(c, traits_type::eof()))
          *this->gptr() = traits_type::to_char_type(c);
        return traits_type::not_eof(c);
      }
      else
         return traits_type::eof();
    }

  template <typename C, typename T>
    std::streamsize
    basic_pstreambuf<C,T>::showmanyc()
    {
      int avail = 0;
      if (sizeof(char_type) == 1)
        avail = fill_buffer(true) ? this->egptr() - this->gptr() : -1;
#ifdef FIONREAD
      else
      {
        if (::ioctl(rpipe(), FIONREAD, &avail) == -1)
          avail = -1;
        else if (avail)
          avail /= sizeof(char_type);
      }
#endif
      return std::streamsize(avail);
    }

  /**
   * @return  true if the buffer was filled, false otherwise.
   */
  template <typename C, typename T>
    bool
    basic_pstreambuf<C,T>::fill_buffer(bool non_blocking)
    {
      const std::streamsize pb1 = this->gptr() - this->eback();
      const std::streamsize pb2 = pbsz;
      const std::streamsize npb = std::min(pb1, pb2);

      char_type* const rbuf = rbuffer();

      if (npb)
        traits_type::move(rbuf + pbsz - npb, this->gptr() - npb, npb);

      std::streamsize rc = -1;

      if (non_blocking)
      {
        const int flags = ::fcntl(rpipe(), F_GETFL);
        if (flags != -1)
        {
          const bool blocking = !(flags & O_NONBLOCK);
          if (blocking)
            ::fcntl(rpipe(), F_SETFL, flags | O_NONBLOCK);  // set non-blocking

          error_ = 0;
          rc = read(rbuf + pbsz, bufsz - pbsz);

          if (rc == -1 && error_ == EAGAIN)  // nothing available
            rc = 0;
          else if (rc == 0)  // EOF
            rc = -1;

          if (blocking)
            ::fcntl(rpipe(), F_SETFL, flags); // restore
        }
      }
      else
        rc = read(rbuf + pbsz, bufsz - pbsz);

      if (rc > 0 || (rc == 0 && non_blocking))
      {
        this->setg( rbuf + pbsz - npb,
                    rbuf + pbsz,
                    rbuf + pbsz + rc );
        return true;
      }
      else
      {
        this->setg(NULL, NULL, NULL);
        return false;
      }
    }

  /**
   * Writes up to @a n characters to the pipe from the buffer @a s.
   *
   * @param   s  character buffer.
   * @param   n  buffer length.
   * @return  the number of characters written.
   */
  template <typename C, typename T>
    inline std::streamsize
    basic_pstreambuf<C,T>::write(const char_type* s, std::streamsize n)
    {
      std::streamsize nwritten = 0;
      if (wpipe() >= 0)
      {
        nwritten = ::write(wpipe(), s, n * sizeof(char_type));
        if (nwritten == -1)
          error_ = errno;
        else
          nwritten /= sizeof(char_type);
      }
      return nwritten;
    }

  /**
   * Reads up to @a n characters from the pipe to the buffer @a s.
   *
   * @param   s  character buffer.
   * @param   n  buffer length.
   * @return  the number of characters read.
   */
  template <typename C, typename T>
    inline std::streamsize
    basic_pstreambuf<C,T>::read(char_type* s, std::streamsize n)
    {
      std::streamsize nread = 0;
      if (rpipe() >= 0)
      {
        nread = ::read(rpipe(), s, n * sizeof(char_type));
        if (nread == -1)
          error_ = errno;
        else
          nread /= sizeof(char_type);
      }
      return nread;
    }

  /** @return a reference to the output file descriptor */
  template <typename C, typename T>
    inline pstreams::fd_type&
    basic_pstreambuf<C,T>::wpipe()
    {
      return wpipe_;
    }

  /** @return a reference to the active input file descriptor */
  template <typename C, typename T>
    inline pstreams::fd_type&
    basic_pstreambuf<C,T>::rpipe()
    {
      return rpipe_[rsrc_];
    }

  /** @return a reference to the specified input file descriptor */
  template <typename C, typename T>
    inline pstreams::fd_type&
    basic_pstreambuf<C,T>::rpipe(buf_read_src which)
    {
      return rpipe_[which];
    }

  /** @return a pointer to the start of the active input buffer area. */
  template <typename C, typename T>
    inline typename basic_pstreambuf<C,T>::char_type*
    basic_pstreambuf<C,T>::rbuffer()
    {
      return rbuffer_[rsrc_];
    }


  /*
   * member definitions for pstream_common
   */

  /**
   * @class pstream_common
   * Abstract Base Class providing common functionality for basic_ipstream,
   * basic_opstream and basic_pstream.
   * pstream_common manages the basic_pstreambuf stream buffer that is used
   * by the derived classes to initialise an iostream class.
   */

  /** Creates an uninitialised stream. */
  template <typename C, typename T>
    inline
    pstream_common<C,T>::pstream_common()
    : std::basic_ios<C,T>(NULL)
    , command_()
    , buf_()
    {
      this->std::basic_ios<C,T>::rdbuf(&buf_);
    }

  /**
   * Initialises the stream buffer by calling
   * do_open( @a command , @a mode )
   *
   * @param cmd   a string containing a shell command.
   * @param mode  the I/O mode to use when opening the pipe.
   * @see   do_open(const std::string&, pmode)
   */
  template <typename C, typename T>
    inline
    pstream_common<C,T>::pstream_common(const std::string& cmd, pmode mode)
    : std::basic_ios<C,T>(NULL)
    , command_(cmd)
    , buf_()
    {
      this->std::basic_ios<C,T>::rdbuf(&buf_);
      do_open(cmd, mode);
    }

  /**
   * Initialises the stream buffer by calling
   * do_open( @a file , @a argv , @a mode )
   *
   * @param file  a string containing the pathname of a program to execute.
   * @param argv  a vector of argument strings passed to the new program.
   * @param mode  the I/O mode to use when opening the pipe.
   * @see do_open(const std::string&, const argv_type&, pmode)
   */
  template <typename C, typename T>
    inline
    pstream_common<C,T>::pstream_common( const std::string& file,
                                         const argv_type& argv,
                                         pmode mode )
    : std::basic_ios<C,T>(NULL)
    , command_(file)
    , buf_()
    {
      this->std::basic_ios<C,T>::rdbuf(&buf_);
      do_open(file, argv, mode);
    }

  /**
   * This is a pure virtual function to make @c pstream_common abstract.
   * Because it is the destructor it will be called by derived classes
   * and so must be defined.  It is also protected, to discourage use of
   * the PStreams classes through pointers or references to the base class.
   *
   * @sa If defining a pure virtual seems odd you should read
   * http://www.gotw.ca/gotw/031.htm (and the rest of the site as well!)
   */
  template <typename C, typename T>
    inline
    pstream_common<C,T>::~pstream_common()
    {
    }

  /**
   * Calls rdbuf()->open( @a command , @a mode )
   * and sets @c failbit on error.
   *
   * @param cmd   a string containing a shell command.
   * @param mode  the I/O mode to use when opening the pipe.
   * @see   basic_pstreambuf::open(const std::string&, pmode)
   */
  template <typename C, typename T>
    inline void
    pstream_common<C,T>::do_open(const std::string& cmd, pmode mode)
    {
      if (!buf_.open((command_=cmd), mode))
        this->setstate(std::ios_base::failbit);
    }

  /**
   * Calls rdbuf()->open( @a file, @a  argv, @a mode )
   * and sets @c failbit on error.
   *
   * @param file  a string containing the pathname of a program to execute.
   * @param argv  a vector of argument strings passed to the new program.
   * @param mode  the I/O mode to use when opening the pipe.
   * @see   basic_pstreambuf::open(const std::string&, const argv_type&, pmode)
   */
  template <typename C, typename T>
    inline void
    pstream_common<C,T>::do_open( const std::string& file,
                                  const argv_type& argv,
                                  pmode mode )
    {
      if (!buf_.open((command_=file), argv, mode))
        this->setstate(std::ios_base::failbit);
    }

  /** Calls rdbuf->close() and sets @c failbit on error. */
  template <typename C, typename T>
    inline void
    pstream_common<C,T>::close()
    {
      if (!buf_.close())
        this->setstate(std::ios_base::failbit);
    }

  /**
   * @return  rdbuf()->is_open().
   * @see     basic_pstreambuf::is_open()
   */
  template <typename C, typename T>
    inline bool
    pstream_common<C,T>::is_open() const
    {
      return buf_.is_open();
    }

  /** @return a pointer to the private stream buffer member. */
  // TODO  document behaviour if buffer replaced.
  template <typename C, typename T>
    inline typename pstream_common<C,T>::streambuf_type*
    pstream_common<C,T>::rdbuf() const
    {
      return const_cast<streambuf_type*>(&buf_);
    }


#if REDI_EVISCERATE_PSTREAMS
  /**
   * @def REDI_EVISCERATE_PSTREAMS
   * If this macro has a non-zero value then certain internals of the
   * @c basic_pstreambuf template class are exposed. In general this is
   * a Bad Thing, as the internal implementation is largely undocumented
   * and may be subject to change at any time, so this feature is only
   * provided because it might make PStreams useful in situations where
   * it is necessary to do Bad Things.
   */

  /**
   * @warning  This function exposes the internals of the stream buffer and
   *           should be used with caution. It is the caller's responsibility
   *           to flush streams etc. in order to clear any buffered data.
   *           The POSIX.1 function <b>fdopen</b>(3) is used to obtain the
   *           @c FILE pointers from the streambuf's private file descriptor
   *           members so consult your system's documentation for
   *           <b>fdopen</b>(3).
   *
   * @param   in    A FILE* that will refer to the process' stdin.
   * @param   out   A FILE* that will refer to the process' stdout.
   * @param   err   A FILE* that will refer to the process' stderr.
   * @return  An OR of zero or more of @c pstdin, @c pstdout, @c pstderr.
   *
   * For each open stream shared with the child process a @c FILE* is
   * obtained and assigned to the corresponding parameter. For closed
   * streams @c NULL is assigned to the parameter.
   * The return value can be tested to see which parameters should be
   * @c !NULL by masking with the corresponding @c pmode value.
   *
   * @see <b>fdopen</b>(3)
   */
  template <typename C, typename T>
    std::size_t
    basic_pstreambuf<C,T>::fopen(FILE*& in, FILE*& out, FILE*& err)
    {
      in = out = err = NULL;
      std::size_t open_files = 0;
      if (wpipe() > -1)
      {
        if ((in = ::fdopen(wpipe(), "w")))
        {
            open_files |= pstdin;
        }
      }
      if (rpipe(rsrc_out) > -1)
      {
        if ((out = ::fdopen(rpipe(rsrc_out), "r")))
        {
            open_files |= pstdout;
        }
      }
      if (rpipe(rsrc_err) > -1)
      {
        if ((err = ::fdopen(rpipe(rsrc_err), "r")))
        {
            open_files |= pstderr;
        }
      }
      return open_files;
    }

  /**
   *  @warning This function exposes the internals of the stream buffer and
   *  should be used with caution.
   *
   *  @param  in   A FILE* that will refer to the process' stdin.
   *  @param  out  A FILE* that will refer to the process' stdout.
   *  @param  err  A FILE* that will refer to the process' stderr.
   *  @return A bitwise-or of zero or more of @c pstdin, @c pstdout, @c pstderr.
   *  @see    basic_pstreambuf::fopen()
   */
  template <typename C, typename T>
    inline std::size_t
    pstream_common<C,T>::fopen(FILE*& fin, FILE*& fout, FILE*& ferr)
    {
      return buf_.fopen(fin, fout, ferr);
    }

#endif // REDI_EVISCERATE_PSTREAMS


} // namespace redi


#endif  // REDI_PSTREAM_H_SEEN
#endif  // WIN32

/* pstream END */
/* main START */

inline string formatANISOU(const string &inputString)
{
    string result=Trim(inputString," ");
    size_t found=result.find_first_of('.');
    stringstream buf;
    if (found==string::npos)
    {
        result+=".0000";
        found=result.find_first_of('.');
    }
    string post_decimal=result.substr(found+1);
    int i;
    long int second=atoi(post_decimal.c_str());
    if (result[0]=='-') second=-second;
    if (post_decimal.size()>4)
    {
        double anisou_dbl=atof(result.c_str())*10000;
        long int anisou_int=round(anisou_dbl);
        if (rstrip(post_decimal.substr(4),"0")=="5" &&
            (post_decimal[3]-'0') % 2 ==0)  anisou_int=(long int)(anisou_dbl);
        buf<<right<<setw(7)<<anisou_int<<flush;
        result=buf.str();
        buf.str(string());
        post_decimal.clear();
        return result.substr(0,7);
    }
    else if (post_decimal.size()<4)
        for (i=0;i<4-post_decimal.size();i++) second*=10;
    
    string pre_decimal=result.substr(0,found);
    long int first=atoi(pre_decimal.c_str());
    buf<<right<<setw(7)<<first*10000+second<<flush;
    result=buf.str();
    buf.str(string());
    //if (result.size()>7) result=result.substr(7-result.size());
    if (result.size()>7) result=result.substr(0,7);
    pre_decimal.clear();
    post_decimal.clear();
    return result;
}

int read_semi_colon(vector<string> &line_vec, const int fields, int l,
    const vector<string> &lines, vector<string> &line_append_vec, string &line,
    const bool ignore_quotation=false,const bool add_space=false)
{
    int i;
    if (line_vec.size() && StartsWith(line_vec[0],";"))
    {
        l--;
        clear_line_vec(line_vec);
    }
    while (line_vec.size()<fields)
    {
        l++;
        if (StartsWith(lines[l],";"))
        {
            line="";
            while (l<lines.size())
            {
                if (StartsWith(lines[l],";"))
                {
                    if (Trim(lines[l])==";") break;
                    else line+=lines[l].substr(1);
                }
                else line+=lines[l];
                if (add_space) line+=' ';
                l++;
            }
            line_vec.push_back(line);
        }
        else
        {
            Split(lines[l],line_append_vec,' ',ignore_quotation);
            for (i=0;i<line_append_vec.size();i++)
                line_vec.push_back(line_append_vec[i]);
            clear_line_vec(line_append_vec);
        }
    }
    return l;
}

string BeEM(string &infile, string &pdbid, const int read_seqres,
    const int read_dbref, const int do_gzip, const int do_upper,
    const long int maxatom, const int outfmt, const string &idmap,
    const vector<string>&ccd3_vec, const vector<string>&outputChain_vec)
{
    stringstream buf;
    vector<string> lines;
    Split(infile,lines,'\n',true); 
    buf.str(string());
    if (lines.size()<=1)
    {
        cerr<<"ERROR! Empty structure "<<infile<<endl;
        vector<string>().swap(lines);
        return 0;
    }

    /* parse PDB ID
     * HEADER, AUTHOR, JRNL, CRYST1, SCALEn */
    vector<string> ccd5_vec;
    map<string,string> ccd5_map;
    string pdbx_keywords="";
    string recvd_initial_deposition_date="";
    string revision_date="";
    bool loop_=false;

    map<string,int> _audit_author;
    map<string,int> _citation_author;
    map<string,int> _citation;
    map<string,int> _cell;
    map<string,int> fract_transf_;
    map<string,int> _atom_site;
    map<string,int> _symmetry;
    map<string,int> _struct_keywords;
    map<string,int> _pdbx_database_status;
    map<string,int> _pdbx_audit_revision_history;
    map<string,int> _entity_poly_seq;
    map<string,int> _entity_poly;
    map<string,int> _struct_ref;
    map<string,int> _struct_ref_seq;
    string _citation_title="";
    string _citation_pdbx_database_id_PubMed="";
    string _citation_pdbx_database_id_DOI="";
    string _citation_journal_abbrev="";
    string _citation_journal_volume="";
    string _citation_page_first="";
    string _citation_year="";
    string _citation_journal_id_ASTM="";
    string _citation_country="";
    string _citation_journal_id_ISSN="";

    string entity_id="";
    string mon_id="";
    string pdbx_strand_id="";

    string group_PDB  ="ATOM"; // (ATOM/HETATM)
    string type_symbol="C";    // (element symbol)
    string atom_id    ="CA";   // auth_atom_id, label_atom_id (atom name)
    string alt_id     =" ";    // auth_alt_id, label_alt_id
                               // (alternative location indicator)
    string comp_id    ="UNK";  // auth_comp_id, label_comp_id (residue name)
    string asym_id    ="A"; // auth_asym_id, label_asym_id (chain ID)
    string seq_id     ="   1"; // label_seq_id, auth_seq_id (residue index)
    string pdbx_PDB_ins_code=" ";// (insertion code)
    string Cartn_x    ="   0.000"; 
    string Cartn_y    ="   0.000"; 
    string Cartn_z    ="   0.000"; 
    string occupancy  ="  1.00";
    string B_iso_or_equiv="  0.00";   // Bfactor
    string pdbx_formal_charge="  ";
    string pdbx_PDB_model_num="   1"; // model index
    vector <string> model_num_vec(1,pdbx_PDB_model_num);
    string U11="  10000";
    string U12="      0";
    string U13="      0";
    string U22="  10000";
    string U23="      0";
    string U33="  10000";

    map<string,string> anisou_map;
    map<string,size_t> chainAtomNum_map;
    map<string,size_t> chainHydrNum_map;
    vector<string> chainID_vec;
    vector<pair<string,string> > atomLine_vec;
    vector<pair<string,string> > ligLine_vec;
    vector<pair<string,string> > hohLine_vec;
    vector<string> seqres_vec;
    map<int,vector<string> > seqres_mat;
    vector<string> entity2strand;
    
    map<string,string> accession2db_name;
    map<string,string> accession2db_code;
    vector<vector<string> > dbref_mat;
    vector<string> dbref_vec(13,""); // pdbx_PDB_id_code, pdbx_strand_id,
    // seq_align_beg, pdbx_seq_align_beg_ins_code
    // seq_align_beg, pdbx_seq_align_end_ins_code
    // db_name, pdbx_db_accession, db_code
    // db_align_beg, pdbx_db_align_beg_ins_code
    // db_align_end, pdbx_db_align_end_ins_code
    dbref_vec[3]=dbref_vec[5]=dbref_vec[10]=dbref_vec[12]=" ";
    string db_code;
    string db_name;
    //string pdbx_PDB_id_code;
    string pdbx_db_accession;
    //string seq_align_beg;
    //string seq_align_end;
    //string pdbx_seq_align_beg_ins_code=" ";
    //string pdbx_seq_align_end_ins_code=" ";
    //string db_align_beg;
    //string db_align_end;
    //string pdbx_db_align_beg_ins_code=" ";
    //string pdbx_db_align_end_ins_code=" ";

    size_t l;
    int i,j;
    string line;
    vector<string> line_vec;
    vector<string> line_append_vec;
    vector<string> author_vec;
    vector<string> citation_author_vec;
    vector<string> cryst1_vec(8,"");
    vector<string> scale_vec(4,"");
    vector<vector<string> > scale_mat(3,scale_vec);
    for (l=0;l<lines.size();l++)
    {
        line=lines[l];
        
        if (_atom_site.size() && !StartsWith(line,"_atom_site"))
             Split(line,line_vec,' ',true);
        else Split(line,line_vec,' ');
        //cout<<"["<<l<<"] "<<line<<endl;
        if (line_vec.size()==0) continue;
        else if (line_vec.size()==1 && line_vec[0]=="#")
        {
            if (entity_id.size())
            {
                if (_entity_poly_seq.size())
                {
                    seqres_mat[atoi(entity_id.c_str())]=seqres_vec;
                    for (i=0;i<seqres_vec.size();i++) seqres_vec[i].clear();
                    seqres_vec.clear();
                }
                else if (read_seqres && pdbx_strand_id.size())
                {
                    i=atoi(entity_id.c_str());
                    while (entity2strand.size()<=i) entity2strand.push_back("");
                    entity2strand[i]=pdbx_strand_id;
                }
                entity_id.clear();
            }
            pdbx_strand_id.clear();
            
            if (pdbx_db_accession.size())
            {
                if (db_name.size()) accession2db_name[pdbx_db_accession]=db_name;
                if (db_code.size()) accession2db_code[pdbx_db_accession]=db_code;
                pdbx_db_accession.clear();
            }
            db_name.clear();
            db_code.clear();

            if (read_dbref && dbref_vec[1].size() && dbref_vec[7].size())
            {
                if (dbref_vec[0].size()==0) dbref_vec[0]=Upper(pdbid);
                for (i=0;i<dbref_vec.size();i++)
                    if (dbref_vec[i]=="?" || dbref_vec[i]==".") dbref_vec[i]=" ";
                dbref_mat.push_back(dbref_vec);
                for (i=0;i<dbref_vec.size();i++) dbref_vec[i]="";
                dbref_vec[3]=dbref_vec[5]=dbref_vec[10]=dbref_vec[12]=" ";
            }

            _audit_author.clear();
            _citation_author.clear();
            _citation.clear();
            _cell.clear();
            fract_transf_.clear();
            _atom_site.clear();
            _symmetry.clear();
            _struct_keywords.clear();
            _pdbx_database_status.clear();
            _pdbx_audit_revision_history.clear();
            _entity_poly_seq.clear();
            _entity_poly.clear();
            _struct_ref.clear();
            _struct_ref_seq.clear();
            loop_=false;
        }
        else if (line_vec.size()==1 && line_vec[0]=="loop_")
            loop_=true;
        else if (pdbid.size()==0 && l==0 && StartsWith(line,"data_"))
            pdbid=Lower(line.substr(5));
        else if (pdbid.size()==0 && line_vec.size()>1 && line_vec[0]=="_entry.id")
            pdbid=Lower(line_vec[1]);
        else if (StartsWith(line,"_struct_keywords"))
        {
            if (loop_)
            {
                j=_struct_keywords.size();
                line=line_vec[0];
                _struct_keywords[line]=j;
            }
            else if (line_vec.size()>1)
            {
                if (line_vec[0]=="_struct_keywords.pdbx_keywords")
                    pdbx_keywords+=Trim(line_vec[1],"'\"");
            }
        }
        else if (_struct_keywords.size() && 
            _struct_keywords.count("_struct_keywords.pdbx_keywords"))
        {
            pdbx_keywords+=Trim(line_vec[_struct_keywords[
                "_struct_keywords.pdbx_keywords"]],"'\"");
        }
        else if (StartsWith(line,"_pdbx_database_status."))
        {
            if (loop_)
            {
                j=_pdbx_database_status.size();
                line=line_vec[0];
                clear_line_vec(line_vec);
                Split(line,line_vec,'.');
                if (line_vec.size()>1)
                {
                    line=line_vec[1];
                    _pdbx_database_status[line]=j;
                }

            }
            else if (line_vec.size()>1)
            {
                if (line_vec[0]=="_pdbx_database_status.recvd_initial_deposition_date")
                    recvd_initial_deposition_date=line_vec[1];
            }
        }
        else if (_pdbx_database_status.size() && 
            _pdbx_database_status.count("recvd_initial_deposition_date"))
        {
            recvd_initial_deposition_date=line_vec[
                _pdbx_database_status["recvd_initial_deposition_date"]];
        }
        else if (read_dbref && StartsWith(line,"_struct_ref."))
        {
            if (loop_)
            {
                j=_struct_ref.size();
                line=line_vec[0];
                clear_line_vec(line_vec);
                Split(line,line_vec,'.');
                if (line_vec.size()>1)
                {
                    line=line_vec[1];
                    _struct_ref[line]=j;
                }
            }
            else
            {
                l=read_semi_colon(line_vec, 2, l, lines, line_append_vec, line);
                if (line_vec[0]=="_struct_ref.db_name")
                    db_name=line_vec[1];
                else if (line_vec[0]=="_struct_ref.db_code")
                    db_code=line_vec[1];
                else if (line_vec[0]=="_struct_ref.pdbx_db_accession")
                    pdbx_db_accession=line_vec[1];
            }
        }
        else if (read_dbref && _struct_ref.count("pdbx_db_accession"))
        {
            l=read_semi_colon(line_vec, _struct_ref.size(),
                l, lines, line_append_vec, line);
            pdbx_db_accession=line_vec[_struct_ref["pdbx_db_accession"]];
            if (_struct_ref.count("db_name")) accession2db_name[
                pdbx_db_accession]=line_vec[_struct_ref["db_name"]];
            if (_struct_ref.count("db_code")) accession2db_code[
                pdbx_db_accession]=line_vec[_struct_ref["db_code"]];
        }
        else if (read_dbref && StartsWith(line,"_struct_ref_seq."))
        {
            if (loop_)
            {
                j=_struct_ref_seq.size();
                line=line_vec[0];
                clear_line_vec(line_vec);
                Split(line,line_vec,'.');
                if (line_vec.size()>1)
                {
                    line=line_vec[1];
                    _struct_ref_seq[line]=j;
                }
            }
            else
            {
                l=read_semi_colon(line_vec, 2, l, lines, line_append_vec, line);
                if (line_vec[0]=="_struct_ref_seq.pdbx_PDB_id_code")
                    dbref_vec[0]=line_vec[1];
                else if (line_vec[0]=="_struct_ref_seq.pdbx_strand_id")
                    dbref_vec[1]=line_vec[1];
                else if (line_vec[0]=="_struct_ref_seq.seq_align_beg" && 
                    dbref_vec[2].size()==0) dbref_vec[2]=line_vec[1];
                else if (line_vec[0]=="_struct_ref_seq.pdbx_auth_seq_align_beg")
                    dbref_vec[2]=line_vec[1];
                else if (line_vec[0]=="_struct_ref_seq.pdbx_seq_align_beg_ins_code")
                    dbref_vec[3]=line_vec[1];
                else if (line_vec[0]=="_struct_ref_seq.seq_align_end" && 
                    dbref_vec[4].size()==0) dbref_vec[4]=line_vec[1];
                else if (line_vec[0]=="_struct_ref_seq.pdbx_auth_seq_align_end")
                    dbref_vec[4]=line_vec[1];
                else if (line_vec[0]=="_struct_ref_seq.pdbx_seq_align_end_ins_code")
                    dbref_vec[5]=line_vec[1];
                else if (line_vec[0]=="_struct_ref_seq.pdbx_db_accession")
                    dbref_vec[7]=line_vec[1];
                else if (line_vec[0]=="_struct_ref_seq.db_align_beg")
                    dbref_vec[9]=line_vec[1];
                else if (line_vec[0]=="_struct_ref_seq.pdbx_db_align_beg_ins_code"
                    && line_vec[1]!="?" && line_vec[1]!=".")
                    dbref_vec[10]=line_vec[1];
                else if (line_vec[0]=="_struct_ref_seq.db_align_end")
                    dbref_vec[11]=line_vec[1];
                else if (line_vec[0]=="_struct_ref_seq.pdbx_db_align_end_ins_code"
                    && line_vec[1]!="?" && line_vec[1]!=".")
                    dbref_vec[12]=line_vec[1];
            }
        }
        
        else if (read_dbref && _struct_ref_seq.count("pdbx_strand_id") &&
            _struct_ref_seq.count("pdbx_db_accession"))
        {
            l=read_semi_colon(line_vec, _struct_ref_seq.size(),
                l, lines, line_append_vec, line);
            if (_struct_ref_seq.count("pdbx_PDB_id_code"))
                dbref_vec[0]=line_vec[_struct_ref_seq["pdbx_PDB_id_code"]];
            else dbref_vec[0]=Upper(pdbid);
            dbref_vec[1]=line_vec[_struct_ref_seq["pdbx_strand_id"]];
            if (_struct_ref_seq.count("pdbx_auth_seq_align_beg"))
                dbref_vec[2]=line_vec[_struct_ref_seq["pdbx_auth_seq_align_beg"]];
            else if (_struct_ref_seq.count("seq_align_beg"))
                dbref_vec[2]=line_vec[_struct_ref_seq["seq_align_beg"]];
            if (_struct_ref_seq.count("pdbx_seq_align_beg_ins_code"))
                dbref_vec[3]=line_vec[_struct_ref_seq["pdbx_seq_align_beg_ins_code"]];
            if (_struct_ref_seq.count("pdbx_auth_seq_align_end"))
                dbref_vec[4]=line_vec[_struct_ref_seq["pdbx_auth_seq_align_end"]];
            else if (_struct_ref_seq.count("seq_align_end"))
                dbref_vec[4]=line_vec[_struct_ref_seq["seq_align_end"]];
            if (_struct_ref_seq.count("pdbx_seq_align_end_ins_code"))
                dbref_vec[5]=line_vec[_struct_ref_seq["pdbx_seq_align_end_ins_code"]];
            dbref_vec[7]=line_vec[_struct_ref_seq["pdbx_db_accession"]];
            if (_struct_ref_seq.count("db_align_beg"))
                dbref_vec[9]=line_vec[_struct_ref_seq["db_align_beg"]];
            if (_struct_ref_seq.count("pdbx_db_align_beg_ins_code"))
                dbref_vec[10]=line_vec[_struct_ref_seq["pdbx_db_align_beg_ins_code"]];
            if (_struct_ref_seq.count("db_align_end"))
                dbref_vec[11]=line_vec[_struct_ref_seq["db_align_end"]];
            if (_struct_ref_seq.count("pdbx_db_align_end_ins_code"))
                dbref_vec[12]=line_vec[_struct_ref_seq["pdbx_db_align_end_ins_code"]];
            
            for (i=0;i<dbref_vec.size();i++)
                if (dbref_vec[i]=="?" || dbref_vec[i]==".") dbref_vec[i]=" ";
            dbref_mat.push_back(dbref_vec);
            for (i=0;i<dbref_vec.size();i++) dbref_vec[i]="";
            dbref_vec[3]=dbref_vec[5]=dbref_vec[10]=dbref_vec[12]=" ";
        }
        else if (read_seqres && StartsWith(line,"_entity_poly."))
        {
            if (loop_)
            {
                j=_entity_poly.size();
                line=line_vec[0];
                clear_line_vec(line_vec);
                Split(line,line_vec,'.');
                if (line_vec.size()>1)
                {
                    line=line_vec[1];
                    _entity_poly[line]=j;
                }
            }
            else
            {
                l=read_semi_colon(line_vec, 2, l, lines, line_append_vec, line);
                if (line_vec[0]=="_entity_poly.entity_id")
                    entity_id=line_vec[1];
                else if (line_vec[0]=="_entity_poly.pdbx_strand_id")
                    pdbx_strand_id=line_vec[1];
            }
        }
        else if (read_seqres && _entity_poly.count("entity_id") && 
            _entity_poly.count("pdbx_strand_id"))
        {
            l=read_semi_colon(line_vec, _entity_poly.size(),
                l, lines, line_append_vec, line);
            entity_id=line_vec[_entity_poly["entity_id"]];
            i=atoi(entity_id.c_str());
            while (entity2strand.size()<=i) entity2strand.push_back("");
            entity2strand[i]=line_vec[_entity_poly["pdbx_strand_id"]];
        }
        else if (read_seqres && StartsWith(line,"_entity_poly_seq.") && loop_)
        {
            j=_entity_poly_seq.size();
            line=line_vec[0];
            clear_line_vec(line_vec);
            Split(line,line_vec,'.');
            if (line_vec.size()>1)
            {
                line=line_vec[1];
                _entity_poly_seq[line]=j;
            }
        }
        else if (read_seqres && _entity_poly_seq.size() && 
            _entity_poly_seq.count("entity_id") && _entity_poly_seq.count("mon_id"))
        {
            mon_id   =line_vec[_entity_poly_seq["mon_id"]];
            if      (mon_id.size()==1) mon_id="  "+mon_id;
            else if (mon_id.size()==2) mon_id=" "+mon_id;
            else if (mon_id.size()>3)  mon_id=mon_id.substr(0,3);
            if (entity_id!=line_vec[_entity_poly_seq["entity_id"]])
            {
                if (entity_id.size())
                {
                    seqres_mat[atoi(entity_id.c_str())]=seqres_vec;
                    for (i=0;i<seqres_vec.size();i++) seqres_vec[i].clear();
                    seqres_vec.clear();
                }
                entity_id=line_vec[_entity_poly_seq["entity_id"]];
            }
            seqres_vec.push_back(mon_id);
        }
        else if (StartsWith(line,"_pdbx_audit_revision_history."))
        {
            if (loop_)
            {
                j=_pdbx_audit_revision_history.size();
                line=line_vec[0];
                clear_line_vec(line_vec);
                Split(line,line_vec,'.');
                if (line_vec.size()>1)
                {
                    line=line_vec[1];
                    _pdbx_audit_revision_history[line]=j;
                }

            }
            else if (line_vec.size()>1)
            {
                if (line_vec[0]=="_pdbx_audit_revision_history.revision_date")
                    revision_date=line_vec[1];
            }
        }
        else if (_pdbx_audit_revision_history.size() && 
            revision_date.size()==0 &&
            _pdbx_audit_revision_history.count("revision_date"))
        {
            revision_date=line_vec[_pdbx_audit_revision_history["revision_date"]];
        }
        else if (StartsWith(line,"_citation."))
        {
            if (loop_)
            {
                line=line_vec[0];
                clear_line_vec(line_vec);
                Split(line,line_vec,'.');
                if (line_vec.size()>1)
                {
                    j=_citation.size();
                    line=line_vec[1];
                    _citation[line]=j;
                }
            }
            else
            {
                l=read_semi_colon(line_vec, 2, l, lines, line_append_vec, line,false,true);
                if      (line_vec[0]=="_citation.title")
                    _citation_title=Trim(line_vec[1],"'\"");
                else if (line_vec[0]=="_citation.pdbx_database_id_PubMed")
                    _citation_pdbx_database_id_PubMed=line_vec[1];
                else if (line_vec[0]=="_citation.pdbx_database_id_DOI")
                    _citation_pdbx_database_id_DOI=line_vec[1];
                else if (line_vec[0]=="_citation.journal_abbrev")
                    _citation_journal_abbrev=Trim(line_vec[1],"'\"");
                else if (line_vec[0]=="_citation.journal_volume")
                    _citation_journal_volume=line_vec[1];
                else if (line_vec[0]=="_citation.page_first")
                    _citation_page_first=line_vec[1];
                else if (line_vec[0]=="_citation.year")
                    _citation_year=line_vec[1];
                else if (line_vec[0]=="_citation.journal_id_ASTM")
                    _citation_journal_id_ASTM=Trim(line_vec[1],"'\"");
                else if (line_vec[0]=="_citation.country")
                    _citation_country=Trim(line_vec[1],"'\"");
                else if (line_vec[0]=="_citation.journal_id_ISSN")
                    _citation_journal_id_ISSN=line_vec[1];
            }
        }
        else if (_citation.size() && (_citation.count("id")==0 ||
                 line_vec[_citation["id"]]=="primary"))
        {
            l=read_semi_colon(line_vec, _citation.size(),
                l, lines, line_append_vec, line,false,true);
            if (_citation.count("title"))
                _citation_title=Trim(line_vec[_citation["title"]],"'\"");
            if (_citation.count("pdbx_database_id_PubMed"))
                _citation_pdbx_database_id_PubMed=line_vec[
                _citation["pdbx_database_id_PubMed"]];
            if (_citation.count("pdbx_database_id_DOI"))
                _citation_pdbx_database_id_DOI=line_vec[
                _citation["pdbx_database_id_DOI"]];
            if (_citation.count("journal_abbrev"))
                _citation_journal_abbrev=Trim(line_vec[
                _citation["journal_abbrev"]],"'\"");
            if (_citation.count("journal_volume"))
                _citation_journal_volume=line_vec[_citation["journal_volume"]];
            if (_citation.count("page_first"))
                _citation_page_first=line_vec[_citation["page_first"]];
            if (_citation.count("year"))
                _citation_year=line_vec[_citation["year"]];
            if (_citation.count("journal_id_ASTM"))
                _citation_journal_id_ASTM=Trim(line_vec[
                _citation["journal_id_ASTM"]],"'\"");
            if (_citation.count("country"))
                _citation_country=Trim(line_vec[_citation["country"]],"'\"");
            if (_citation.count("journal_id_ISSN"))
                _citation_journal_id_ISSN=line_vec[_citation["journal_id_ISSN"]];
        }
        else if (StartsWith(line,"_cell."))
        {
            if (loop_)
            {
                j=_cell.size();
                line=line_vec[0];
                _cell[line]=j;
            }
            else if (line_vec.size()>1)
            {
                if      (line_vec[0]=="_cell.length_a")
                    cryst1_vec[0]=formatString(line_vec[1],9,3);
                else if (line_vec[0]=="_cell.length_b")
                    cryst1_vec[1]=formatString(line_vec[1],9,3);
                else if (line_vec[0]=="_cell.length_c")
                    cryst1_vec[2]=formatString(line_vec[1],9,3);
                else if (line_vec[0]=="_cell.angle_alpha")
                    cryst1_vec[3]=formatString(line_vec[1],7,2);
                else if (line_vec[0]=="_cell.angle_beta")
                    cryst1_vec[4]=formatString(line_vec[1],7,2);
                else if (line_vec[0]=="_cell.angle_gamma")
                    cryst1_vec[5]=formatString(line_vec[1],7,2);
                else if (line_vec[0]=="_cell.Z_PDB")
                    cryst1_vec[7]=line_vec[1];
            }
        }
        else if (_cell.size())
        {
            if (_cell.count("_cell.length_a"))    cryst1_vec[0]=
                formatString(line_vec[_cell["_cell.length_a"]],9,3);
            if (_cell.count("_cell.length_b"))    cryst1_vec[1]=
                formatString(line_vec[_cell["_cell.length_b"]],9,3);
            if (_cell.count("_cell.length_c"))    cryst1_vec[2]=
                formatString(line_vec[_cell["_cell.length_c"]],9,3);
            if (_cell.count("_cell.angle_alpha")) cryst1_vec[3]=
                formatString(line_vec[_cell["_cell.angle_alpha"]],7,2);
            if (_cell.count("_cell.angle_beta"))  cryst1_vec[4]=
                formatString(line_vec[_cell["_cell.angle_beta"]],7,2);
            if (_cell.count("_cell.angle_gamma")) cryst1_vec[5]=
                formatString(line_vec[_cell["_cell.angle_gamma"]],7,2);
            if (_cell.count("_cell.Z_PDB"))       cryst1_vec[7]=
                line_vec[_cell["_cell.Z_PDB"]];
        }
        else if (StartsWith(line,"_symmetry"))
        {
            if (loop_)
            {
                j=_symmetry.size();
                line=line_vec[0];
                _symmetry[line]=j;
            }
            else if (line_vec.size()>1)
            {
                if (line_vec[0]=="_symmetry.space_group_name_H-M")
                    cryst1_vec[6]=Trim(line_vec[1],"'\"").substr(0,11);
            }
        }
        else if (_symmetry.size() && 
             _symmetry.count("_symmetry.space_group_name_H-M"))
        {
            cryst1_vec[6]=Trim(line_vec[_symmetry[
                "_symmetry.space_group_name_H-M"]],"'\"").substr(0,11);
        }
        else if (StartsWith(line,"_atom_sites.fract_transf_"))
        {
            if (loop_)
            {
                line=line_vec[0];
                clear_line_vec(line_vec);
                Split(line,line_vec,'.');
                if (line_vec.size()>1)
                {
                    j=fract_transf_.size();
                    line=line_vec[1];
                    fract_transf_[line]=j;
                }
            }
            else if (line_vec.size()>1)
            {
                if      (line_vec[0]=="_atom_sites.fract_transf_matrix[1][1]")
                    scale_mat[0][0]=formatString(line_vec[1],10,6);
                else if (line_vec[0]=="_atom_sites.fract_transf_matrix[1][2]")
                    scale_mat[0][1]=formatString(line_vec[1],10,6);
                else if (line_vec[0]=="_atom_sites.fract_transf_matrix[1][3]")
                    scale_mat[0][2]=formatString(line_vec[1],10,6);
                else if (line_vec[0]=="_atom_sites.fract_transf_matrix[2][1]")
                    scale_mat[1][0]=formatString(line_vec[1],10,6);
                else if (line_vec[0]=="_atom_sites.fract_transf_matrix[2][2]")
                    scale_mat[1][1]=formatString(line_vec[1],10,6);
                else if (line_vec[0]=="_atom_sites.fract_transf_matrix[2][3]")
                    scale_mat[1][2]=formatString(line_vec[1],10,6);
                else if (line_vec[0]=="_atom_sites.fract_transf_matrix[3][1]")
                    scale_mat[2][0]=formatString(line_vec[1],10,6);
                else if (line_vec[0]=="_atom_sites.fract_transf_matrix[3][2]")
                    scale_mat[2][1]=formatString(line_vec[1],10,6);
                else if (line_vec[0]=="_atom_sites.fract_transf_matrix[3][3]")
                    scale_mat[2][2]=formatString(line_vec[1],10,6);
                else if (line_vec[0]=="_atom_sites.fract_transf_vector[1]")
                    scale_mat[0][3]=formatString(line_vec[1],10,5);
                else if (line_vec[0]=="_atom_sites.fract_transf_vector[2]")
                    scale_mat[1][3]=formatString(line_vec[1],10,5);
                else if (line_vec[0]=="_atom_sites.fract_transf_vector[3]")
                    scale_mat[2][3]=formatString(line_vec[1],10,5);
            }
        }
        else if (fract_transf_.size())
        {
            if (fract_transf_.count("fract_transf_matrix[1][1]"))
                scale_mat[0][0]=formatString(line_vec[fract_transf_[
                     "fract_transf_matrix[1][1]"]],10,6);
            if (fract_transf_.count("fract_transf_matrix[1][2]"))
                scale_mat[0][1]=formatString(line_vec[fract_transf_[
                    "fract_transf_matrix[1][2]"]],10,6);
            if (fract_transf_.count("fract_transf_matrix[1][3]"))
                scale_mat[0][2]=formatString(line_vec[fract_transf_[
                    "fract_transf_matrix[1][3]"]],10,6);
            if (fract_transf_.count("fract_transf_matrix[2][1]"))
                scale_mat[1][0]=formatString(line_vec[fract_transf_[
                    "fract_transf_matrix[2][1]"]],10,6);
            if (fract_transf_.count("fract_transf_matrix[2][2]"))
                scale_mat[1][1]=formatString(line_vec[fract_transf_[
                    "fract_transf_matrix[2][2]"]],10,6);
            if (fract_transf_.count("fract_transf_matrix[2][3]"))
                scale_mat[1][2]=formatString(line_vec[fract_transf_[
                    "fract_transf_matrix[2][3]"]],10,6);
            if (fract_transf_.count("fract_transf_matrix[3][1]"))
                scale_mat[2][0]=formatString(line_vec[fract_transf_[
                    "fract_transf_matrix[3][1]"]],10,6);
            if (fract_transf_.count("fract_transf_matrix[3][2]"))
                scale_mat[2][1]=formatString(line_vec[fract_transf_[
                    "fract_transf_matrix[3][2]"]],10,6);
            if (fract_transf_.count("fract_transf_matrix[3][3]"))
                scale_mat[2][2]=formatString(line_vec[fract_transf_[
                    "fract_transf_matrix[3][3]"]],10,6);
            if (fract_transf_.count("fract_transf_vector[1]"))
                scale_mat[0][3]=formatString(line_vec[fract_transf_[
                    "fract_transf_vector[1]"]],10,5);
            if (fract_transf_.count("fract_transf_vector[2]"))
                scale_mat[1][3]=formatString(line_vec[fract_transf_[
                    "fract_transf_vector[2]"]],10,5);
            if (fract_transf_.count("fract_transf_vector[3]"))
                scale_mat[2][3]=formatString(line_vec[fract_transf_[
                    "fract_transf_vector[3]"]],10,5);
        }
        else if (StartsWith(line,"_audit_author."))
        {
            if (loop_)
            {
                line=line_vec[0];
                j=_audit_author.size();
                _audit_author[line]=j;
            }
            else if (StartsWith(line,"_audit_author.name"))
            {
                l=read_semi_colon(line_vec, 2,
                    l, lines, line_append_vec, line);
                line=line_vec[1];
                line=Trim(line,"'\"");
                clear_line_vec(line_vec);
                Split(line,line_vec,',',true);
                if (line_vec.size()>=2) line=lstrip(line_vec[1])+line_vec[0];
                author_vec.push_back(Upper(line));
            }
        }
        else if (_audit_author.size() && 
                 _audit_author.count("_audit_author.name"))
        {
            l=read_semi_colon(line_vec, _audit_author.size(),
                l, lines, line_append_vec, line);
            line=line_vec[_audit_author["_audit_author.name"]];
            line=Trim(line,"'\"");
            clear_line_vec(line_vec);
            Split(line,line_vec,',',true);
            if (line_vec.size()>=2) line=lstrip(line_vec[1])+line_vec[0];
            author_vec.push_back(Upper(line));
        }
        else if (StartsWith(line,"_citation_author."))
        {
            if (loop_)
            {
                line=line_vec[0];
                clear_line_vec(line_vec);
                Split(line,line_vec,'.');
                line=line_vec[1];
                j=_citation_author.size();
                _citation_author[line]=j;
            }
            else if (StartsWith(line,"_citation_author.name") && line_vec.size()>1)
            {
                line=line_vec[1];
                line=Trim(line,"'\"");
                clear_line_vec(line_vec);
                Split(line,line_vec,',',true);
                if (line_vec.size()>=2)
                    line=lstrip(line_vec[1])+line_vec[0];
                citation_author_vec.push_back(line);
            }
        }
        else if (_citation_author.size() && _citation_author.count("name") 
             && (_citation_author.count("citation_id")==0 || 
             line_vec[_citation_author["citation_id"]]=="primary"))
        {
            l=read_semi_colon(line_vec, _citation_author.size(),
                l, lines, line_append_vec, line);
            line=line_vec[_citation_author["name"]];
            line=Trim(line,"'\"");
            clear_line_vec(line_vec);
            Split(line,line_vec,',',true);
            if (line_vec.size()>=2)
                line=lstrip(line_vec[1])+line_vec[0];
            citation_author_vec.push_back(line);
        }
        else if (StartsWith(line,"_atom_site.") || 
                 StartsWith(line,"_atom_site_anisotrop."))
        {
            line=line_vec[0];
            clear_line_vec(line_vec);
            Split(line,line_vec,'.');
            if (line_vec.size()>1)
            {
                line=line_vec[1];
                j=_atom_site.size();
                _atom_site[line]=j;
            }
        }
        else if (_atom_site.size())
        {
            if (_atom_site.count("group_PDB"))
                group_PDB=line_vec[_atom_site["group_PDB"]];
            if (group_PDB=="ATOM") group_PDB="ATOM  ";
            
            if (_atom_site.count("auth_atom_id"))
            {
                atom_id=line_vec[_atom_site["auth_atom_id"]];
                if (atom_id.size()>4 && _atom_site.count("label_atom_id"))
                    atom_id=line_vec[_atom_site["label_atom_id"]];
            }
            else if (_atom_site.count("label_atom_id"))
                atom_id=line_vec[_atom_site["label_atom_id"]];
            else if (_atom_site.count("pdbx_auth_atom_id"))
            {
                atom_id=line_vec[_atom_site["pdbx_auth_atom_id"]];
                if (atom_id.size()>4 && _atom_site.count("pdbx_label_atom_id"))
                    atom_id=line_vec[_atom_site["pdbx_label_atom_id"]];
            }
            else if (_atom_site.count("pdbx_label_atom_id"))
                atom_id=line_vec[_atom_site["pdbx_label_atom_id"]];
            if ((atom_id[0]=='"'  && atom_id[atom_id.size()-1]=='"')||
                (atom_id[0]=='\'' && atom_id[atom_id.size()-1]=='\''))
                atom_id=atom_id.substr(1,atom_id.size()-2);
            atom_id=atom_id.substr(0,4);
            
            if (_atom_site.count("type_symbol"))
                type_symbol=line_vec[_atom_site["type_symbol"]].substr(0,2);
            else type_symbol=lstrip(atom_id,"1234567890 ")[0];

            if (type_symbol.size()==2) while (atom_id.size()<4) atom_id+=' ';
            else
            {
                if (atom_id.size()==1) atom_id+=' ';
                if (atom_id.size()==2) atom_id+=' ';
                if (atom_id.size()==3) atom_id=' '+atom_id;
            }
            if (type_symbol.size()==1) type_symbol=' '+type_symbol;

            if (_atom_site.count("auth_alt_id"))
                alt_id=line_vec[_atom_site["auth_alt_id"]];
            else if (_atom_site.count("label_alt_id"))
                alt_id=line_vec[_atom_site["label_alt_id"]];
            else if (_atom_site.count("pdbx_auth_alt_id"))
                alt_id=line_vec[_atom_site["pdbx_auth_alt_id"]];
            else if (_atom_site.count("pdbx_label_alt_id"))
                alt_id=line_vec[_atom_site["pdbx_label_alt_id"]];
            if (alt_id=="." || alt_id=="?") alt_id=" ";
            else alt_id=alt_id[0];

            if (_atom_site.count("auth_comp_id"))
            {
                comp_id=line_vec[_atom_site["auth_comp_id"]];
                if (comp_id.size()>3 && _atom_site.count("label_comp_id"))
                    comp_id=line_vec[_atom_site["label_comp_id"]];
            }
            else if (_atom_site.count("label_comp_id"))
                comp_id=line_vec[_atom_site["label_comp_id"]];
            else if (_atom_site.count("pdbx_auth_comp_id"))
            {
                comp_id=line_vec[_atom_site["pdbx_auth_comp_id"]];
                if (comp_id.size()>3 && _atom_site.count("pdbx_label_comp_id"))
                    comp_id=line_vec[_atom_site["pdbx_label_comp_id"]];
            }
            else if (_atom_site.count("pdbx_label_comp_id"))
                comp_id=line_vec[_atom_site["pdbx_label_comp_id"]];
            if      (comp_id.size()==1) comp_id="  "+comp_id;
            else if (comp_id.size()==2) comp_id=" "+comp_id;
            if (comp_id.size()>3)
            {
                if (ccd3_vec.size()==0) comp_id=comp_id.substr(0,3);
                else
                {
                    if (ccd5_map.count(comp_id)==0)
                    {
                        ccd5_map[comp_id]=ccd3_vec[ccd5_vec.size() % ccd3_vec.size()];
                        ccd5_vec.push_back(comp_id);
                    }
                    comp_id=ccd5_map[comp_id];
                }
            }

            if (_atom_site.count("auth_asym_id"))
                asym_id=line_vec[_atom_site["auth_asym_id"]];
            else if (_atom_site.count("label_asym_id"))
                asym_id=line_vec[_atom_site["label_asym_id"]];
            else if (_atom_site.count("pdbx_auth_asym_id"))
                asym_id=line_vec[_atom_site["pdbx_auth_asym_id"]];
            else if (_atom_site.count("pdbx_label_asym_id"))
                asym_id=line_vec[_atom_site["pdbx_label_asym_id"]];
            if (asym_id=="." || asym_id=="?") asym_id="_";
            if (outputChain_vec.size() && find(outputChain_vec.begin(),
                outputChain_vec.end(), asym_id)==outputChain_vec.end())
            {
                clear_line_vec(line_vec);
                continue;
            }
            if (asym_id=="_") asym_id=" ";

            if (_atom_site.count("auth_seq_id"))
                seq_id=line_vec[_atom_site["auth_seq_id"]];
            else if (_atom_site.count("label_seq_id"))
                seq_id=line_vec[_atom_site["label_seq_id"]];
            else if (_atom_site.count("pdbx_auth_seq_id"))
                seq_id=line_vec[_atom_site["pdbx_auth_seq_id"]];
            else if (_atom_site.count("pdbx_label_seq_id"))
                seq_id=line_vec[_atom_site["pdbx_label_seq_id"]];
            if (seq_id.size()>=2) seq_id=lstrip(seq_id,"0");
            if (seq_id.size()==3) seq_id=" "+seq_id;
            else if (seq_id.size()==2) seq_id="  "+seq_id;
            else if (seq_id.size()==1) seq_id="   "+seq_id;
            //else if (seq_id.size()>4) seq_id=seq_id.substr(seq_id.size()-4);

            if (_atom_site.count("pdbx_PDB_ins_code"))
                pdbx_PDB_ins_code=line_vec[_atom_site["pdbx_PDB_ins_code"]];
            if (pdbx_PDB_ins_code=="." || pdbx_PDB_ins_code=="?")
                pdbx_PDB_ins_code=" ";
            else pdbx_PDB_ins_code=pdbx_PDB_ins_code[0];
            seq_id+=pdbx_PDB_ins_code;
            if (seq_id.size()>5) seq_id=seq_id.substr(0,5);

            if (_atom_site.count("Cartn_z"))
            {
                Cartn_x=formatString(line_vec[_atom_site["Cartn_x"]],8,3);
                Cartn_y=formatString(line_vec[_atom_site["Cartn_y"]],8,3);
                Cartn_z=formatString(line_vec[_atom_site["Cartn_z"]],8,3);
            }

            if (_atom_site.count("aniso_U[3][3]"))
            {
                U11=formatANISOU(line_vec[_atom_site["aniso_U[1][1]"]]);
                U12=formatANISOU(line_vec[_atom_site["aniso_U[1][2]"]]);
                U13=formatANISOU(line_vec[_atom_site["aniso_U[1][3]"]]);
                U22=formatANISOU(line_vec[_atom_site["aniso_U[2][2]"]]);
                U23=formatANISOU(line_vec[_atom_site["aniso_U[2][3]"]]);
                U33=formatANISOU(line_vec[_atom_site["aniso_U[3][3]"]]);
            }
            else if (_atom_site.count("U[3][3]"))
            {
                U11=formatANISOU(line_vec[_atom_site["U[1][1]"]]);
                U12=formatANISOU(line_vec[_atom_site["U[1][2]"]]);
                U13=formatANISOU(line_vec[_atom_site["U[1][3]"]]);
                U22=formatANISOU(line_vec[_atom_site["U[2][2]"]]);
                U23=formatANISOU(line_vec[_atom_site["U[2][3]"]]);
                U33=formatANISOU(line_vec[_atom_site["U[3][3]"]]);
            }


            if (_atom_site.count("occupancy"))
                occupancy=formatString(line_vec[_atom_site["occupancy"]],6,2);

            if (_atom_site.count("B_iso_or_equiv"))
                B_iso_or_equiv=formatString(line_vec[_atom_site["B_iso_or_equiv"]],6,2);

            if (_atom_site.count("pdbx_formal_charge"))
                pdbx_formal_charge=line_vec[_atom_site["pdbx_formal_charge"]];
            if (pdbx_formal_charge=="." || pdbx_formal_charge=="?")
                pdbx_formal_charge="  ";
            else if (pdbx_formal_charge.size()==1)
            {
                if ('1'<=pdbx_formal_charge[0] && pdbx_formal_charge[0]<='9')
                     pdbx_formal_charge+='+';
                else pdbx_formal_charge+=' ';
            }
            else if (pdbx_formal_charge.size()>2)
                pdbx_formal_charge=pdbx_formal_charge.substr(0,2);
            if (pdbx_formal_charge[0]=='+' || pdbx_formal_charge[0]=='-')
                pdbx_formal_charge=pdbx_formal_charge.substr(1)+
                                   pdbx_formal_charge[0];

            if (_atom_site.count("pdbx_PDB_model_num"))
                pdbx_PDB_model_num=line_vec[_atom_site["pdbx_PDB_model_num"]];
            if (pdbx_PDB_model_num=="." || pdbx_PDB_model_num=="?")
                pdbx_PDB_model_num="   1";
            else if (pdbx_PDB_model_num.size()==1) 
                pdbx_PDB_model_num="   "+pdbx_PDB_model_num;
            else if (pdbx_PDB_model_num.size()==2) 
                pdbx_PDB_model_num="  "+pdbx_PDB_model_num;
            else if (pdbx_PDB_model_num.size()==3) 
                pdbx_PDB_model_num=" "+pdbx_PDB_model_num;
            if (pdbx_PDB_model_num!="   1" && find(model_num_vec.begin(),
                model_num_vec.end(), pdbx_PDB_model_num)==model_num_vec.end())
                model_num_vec.push_back(pdbx_PDB_model_num);

            //if (pdbx_PDB_model_num=="   1" && _atom_site.count("Cartn_z"))
            if (_atom_site.count("Cartn_z"))
            {
/*
COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.
*/
                line=group_PDB+' '+pdbx_PDB_model_num+' '+atom_id+alt_id
                    +comp_id+"  "+seq_id+"   "
                    +Cartn_x+Cartn_y+Cartn_z+occupancy+B_iso_or_equiv
                    +"          "+type_symbol+pdbx_formal_charge;
                if (_atom_site.count("label_seq_id") && 
                    line_vec[_atom_site["label_seq_id"]]==".")
                {
                    if (comp_id=="HOH") 
                         hohLine_vec.push_back(make_pair(line,asym_id));
                    else ligLine_vec.push_back(make_pair(line,asym_id));
                }
                else atomLine_vec.push_back(make_pair(line,asym_id));
                if (pdbx_PDB_model_num=="   1")
                {
                    if (chainAtomNum_map.count(asym_id)==0)
                    {
                        chainID_vec.push_back(asym_id);
                        chainAtomNum_map[asym_id]=1;
                        chainHydrNum_map[asym_id]=0;
                    }
                    else chainAtomNum_map[asym_id]++;
                    chainHydrNum_map[asym_id]+=(type_symbol==" H");
                }
            }
            if (pdbx_PDB_model_num=="   1" && (_atom_site.count("U[3][3]")||
                                         _atom_site.count("aniso_U[3][3]")))
            {
                /*
COLUMNS       DATA  TYPE    FIELD          DEFINITION
-----------------------------------------------------------------
 1 - 6        Record name   "ANISOU"
 7 - 11       Integer       serial         Atom serial number.
13 - 16       Atom          name           Atom name.
17            Character     altLoc         Alternate location indicator
18 - 20       Residue name  resName        Residue name.
22            Character     chainID        Chain identifier.
23 - 26       Integer       resSeq         Residue sequence number.
27            AChar         iCode          Insertion code.
29 - 35       Integer       u[0][0]        U(1,1)
36 - 42       Integer       u[1][1]        U(2,2)
43 - 49       Integer       u[2][2]        U(3,3)
50 - 56       Integer       u[0][1]        U(1,2)
57 - 63       Integer       u[0][2]        U(1,3)
64 - 70       Integer       u[1][2]        U(2,3)
77 - 78       LString(2)    element        Element symbol, right-justified.
79 - 80       LString(2)    charge         Charge on the atom.
                 */
                anisou_map[atom_id+alt_id+comp_id+"  "+seq_id+
                    '\t'+asym_id]=U11+U22+U33+U12+U13+U23;
            }
        }

        /* clean up */
        clear_line_vec(line_vec);
        lines[l].clear();
    }
    lines.clear();

    if (pdbid.size()==0)
    {
        cerr<<"ERROR: no PDB ID in "<<infile<<'\n'
            <<"PDB ID can be specified by option -p=xxxx"<<endl;
        return " ";
    }


    string header1;
    string header2;
    //if (revision_date.size()) recvd_initial_deposition_date=revision_date;
    size_t found;
    if (pdbx_keywords.size() || recvd_initial_deposition_date.size())
    {
        buf<<"HEADER    "<<left<<setw(40)<<Upper(pdbx_keywords.substr(0,40))
            <<recvd_initial_deposition_date.substr(0,10)<<"  XXXX              "<<endl;
        header1=buf.str();
        buf.str(string());
    }
    if (author_vec.size())
    {
        string author_txt=Join(", ",author_vec);
        for (i=0;i<author_vec.size();i++) author_vec[i].clear();
        author_vec.clear();
        Split(author_txt, author_vec,' ',true);
        author_txt.clear();
        int Continuation=0; 
        line="";
        for (i=0;i<author_vec.size();i++)
        {
            if (line.size()==0)
            {
                Continuation++;
                if (Continuation==1) line="AUTHOR    ";
                else
                {
                    buf<<"AUTHOR  "<<right<<setw(2)<<Continuation<<" ";
                    line=buf.str();
                    buf.str(string());
                }
                line+=author_vec[i];
            }
            else if (author_vec[i].size()+line.size()>=79)
            {
                found=author_vec[i].find_first_of('-');
                if (found!=string::npos && found+line.size()<=77)
                {
                    line+=' '+author_vec[i].substr(0,found+1);
                    author_vec[i]=author_vec[i].substr(found+1);
                }
                buf<<left<<setw(80)<<line<<endl;
                header1+=buf.str();
                buf.str(string());
                line="";
                i--;
            }
            else line+=' '+author_vec[i];
        }
        if (line.size())
        {
            buf<<left<<setw(80)<<line<<endl;
            header1+=buf.str();
            buf.str(string());
        }
    }
    if (citation_author_vec.size())
    {
        int Continuation=0; 
        line="";
        for (i=0;i<citation_author_vec.size();i++)
        {
            if (line.size()==0)
            {
                Continuation++;
                if (Continuation==1) line="JRNL        AUTH   ";
                else
                {
                    buf<<"JRNL        AUTH"<<right<<setw(2)<<Continuation<<" ";
                    line=buf.str();
                    buf.str(string());
                }
                if (line.size()+citation_author_vec[i].size()<=79)
                    line+=citation_author_vec[i];
                else
                {
                    Split(citation_author_vec[i],line_vec,' ');
                    for (j=0;j<line_vec.size();j++)
                    {
                        if (line.size()+line_vec[j].size()<=79)
                            line+=line_vec[j]+' ';
                        else break;
                    }
                    buf<<left<<setw(80)<<line<<endl;
                    header1+=buf.str();
                    buf.str(string());
                    citation_author_vec[i]=Join(" ",line_vec,j);
                    clear_line_vec(line_vec);
                    line="";
                    i--;
                }
            }
            else if (citation_author_vec[i].size()+line.size()>=
                77+(i+1==citation_author_vec.size()))
            {
                buf<<left<<setw(80)<<line+","<<endl;
                header1+=buf.str();
                buf.str(string());
                line="";
                i--;
            }
            else line+=", "+citation_author_vec[i];
        }
        if (line.size())
        {
            buf<<left<<setw(80)<<line<<endl;
            header1+=buf.str();
            buf.str(string());
        }
    }
    if (_citation_title.size())
    {
        Split(_citation_title,line_vec,' ',true);
        int Continuation=0; 
        line="";
        for (i=0;i<line_vec.size();i++)
        {
            if (line.size()==0)
            {
                Continuation++;
                if (Continuation==1) line="JRNL        TITL   ";
                else
                {
                    buf<<"JRNL        TITL "<<Continuation<<" ";
                    line=buf.str();
                    buf.str(string());
                }
                if (line.size()+line_vec[i].size()<79)
                    line+=line_vec[i];
                else
                {
                    bool add_word=false;
                    for (j=line_vec[i].size()-2;j>0;j--)
                    {
                        if (line.size()+j>=78 || line_vec[i][j]!=')'||
                            line_vec[i][j+1]=='-') continue;
                        line+=line_vec[i].substr(0,j+1);
                        line_vec[i]=line_vec[i].substr(j+1);
                        add_word=true;
                        i--;
                        break;
                    }
                    if (add_word==false) 
                        line+=line_vec[i].substr(0,78-line_vec[i].size());
                }
            }
            else if (line.size()+line_vec[i].size()>=79)
            {
                if (line_vec[i].size()>60)
                {
                    for (j=line_vec[i].size()-2;j>0;j--)
                    {
                        if (line.size()+j>=78 || line_vec[i][j]!=')'||
                            line_vec[i][j+1]=='-') continue;
                        line+=' '+line_vec[i].substr(0,j+1);
                        line_vec[i]=line_vec[i].substr(j+1);
                        break;
                    }
                }
                buf<<left<<setw(80)<<line<<endl;
                header1+=buf.str();
                buf.str(string());
                line="";
                i--;
            }
            else line+=" "+line_vec[i];
        }
        if (line.size())
        {
            buf<<left<<setw(80)<<line<<endl;
            header1+=buf.str();
            buf.str(string());
        }
        clear_line_vec(line_vec);
    }
    if (_citation_title=="?")                   _citation_title="";
    if (_citation_pdbx_database_id_PubMed=="?") _citation_pdbx_database_id_PubMed="";
    if (_citation_pdbx_database_id_DOI=="?")    _citation_pdbx_database_id_DOI="";
    if (_citation_journal_abbrev=="?")          _citation_journal_abbrev="";
    if (_citation_journal_volume=="?")          _citation_journal_volume="";
    if (_citation_page_first=="?")              _citation_page_first="";
    if (_citation_year=="?")                    _citation_year="";
    if (_citation_journal_id_ASTM=="?")         _citation_journal_id_ASTM="";
    if (_citation_country=="?")                 _citation_country="";
    if (_citation_journal_id_ISSN=="?")         _citation_journal_id_ISSN="";
    if (_citation_journal_abbrev.size())
    {
        Split(_citation_journal_abbrev,line_vec,' ');
        int Continuation=0; 
        line="";
        for (i=0;i<=line_vec.size();i++)
        {
            if (i<line_vec.size() && line.size()==0)
            {
                Continuation++;
                if (Continuation==1) line="JRNL        REF    ";
                else
                {
                    buf<<"JRNL        REF  "<<Continuation<<" ";
                    line=buf.str();
                    buf.str(string());
                }
                line+=line_vec[i];
            }
            else if (i==line_vec.size() || line_vec[i].size()+line.size()>=47)
            {
                if (Continuation==1)
                {
                    if (_citation_journal_volume.size())
                    {
                        while (_citation_journal_volume.size()<4)
                            _citation_journal_volume=' '+_citation_journal_volume;
                        _citation_journal_volume="V."+_citation_journal_volume;
                    }
                    if (_citation_page_first.size()>5) _citation_page_first=
                        _citation_page_first.substr(_citation_page_first.size()-5);
                    buf<<left<<setw(49)<<line
                        <<setw(6)<<left<<_citation_journal_volume.substr(0,6)
                        <<' '<<setw(5)<<right<<Upper(_citation_page_first)
                        <<' '<<left<<setw(18)<<_citation_year<<endl;;
                }
                else buf<<left<<setw(80)<<line<<endl;
                header1+=buf.str();
                buf.str(string());
                line="";
                if (i<line_vec.size()) i--;
            }
            else line+=" "+line_vec[i];
        }
        clear_line_vec(line_vec);
    //}
    //if (_citation_journal_id_ASTM.size() || _citation_country.size() ||
        //_citation_journal_id_ISSN.size())
    //{
        buf<<"JRNL        REFN   "
            <<setw(11)<<right<<_citation_journal_id_ASTM.substr(0,11)
            <<"  "<<setw(7)<<left<<_citation_country.substr(0,7)<<' '
            <<setw(40)<<left<<_citation_journal_id_ISSN<<endl;
        header1+=buf.str();
        buf.str(string());
    }
    if (_citation_pdbx_database_id_PubMed.size())
    {
        buf<<left<<setw(80)<<"JRNL        PMID   "+
            _citation_pdbx_database_id_PubMed<<endl;
        header1+=buf.str();
        buf.str(string());
    }
    if (_citation_pdbx_database_id_DOI.size())
    {
        buf<<left<<setw(80)<<"JRNL        DOI    "+
            Trim(_citation_pdbx_database_id_DOI,"'\"")<<endl;
        header1+=buf.str();
        buf.str(string());
    }
    int cryst1Count=0;
    for (i=0;i<cryst1_vec.size();i++)
    {
        cryst1Count+=cryst1_vec[i].size()>0;
        if (i>=6 && (cryst1_vec[i]=="?" || cryst1_vec[i]==".")) cryst1_vec[i]="";
    }
    if (cryst1Count)
    {
        buf<<"CRYST1"<<setw(9)<<cryst1_vec[0]<<setw(9)<<cryst1_vec[1]<<setw(9)
            <<cryst1_vec[2]<<setw(7)<<cryst1_vec[3]<<setw(7)<<cryst1_vec[4]
            <<setw(7)<<cryst1_vec[5]<<' '<<setw(11)<<left<<cryst1_vec[6]
            <<setw(4)<<right<<cryst1_vec[7]<<"          "<<endl;
        header2+=buf.str();
        buf.str(string());
    }
    int scaleCount=0;
    for (i=0;i<3;i++) for (j=0;j<4;j++) scaleCount+=scale_mat[i][j].size()>0;
    if (scaleCount==12) header2+=
        "SCALE1    "+scale_mat[0][0]+scale_mat[0][1]+scale_mat[0][2]+
        "     "     +scale_mat[0][3]+"                         \n"+
        "SCALE2    "+scale_mat[1][0]+scale_mat[1][1]+scale_mat[1][2]+
        "     "     +scale_mat[1][3]+"                         \n"+
        "SCALE3    "+scale_mat[2][0]+scale_mat[2][1]+scale_mat[2][2]+
        "     "     +scale_mat[2][3]+"                         \n";

    /* parse extra long chain */
    int atomNum=0;
    int SplitNum;
    string key;
    map<string,map<string,int> > SplitChain_map; // asym_id => (atom key => SplitNum)
    map<string,map<string,int> > SplitChainRes_map; // asym_id => (res key => SplitNum)
    map<string,int> SplitChainNum_map; // asym_id => SplitNum
    string res;
    for (i=0;i<chainID_vec.size();i++)
    {
        asym_id=chainID_vec[i];
        if (maxatom<=0 || outfmt==3 || (maxatom>1 && chainAtomNum_map[asym_id]<maxatom))
            continue;
        map<string,int> key_map; // key => SplitNum
        map<string,int>::iterator it;

        atomNum=0;
        SplitNum=0;
        SplitChainNum_map[asym_id]=SplitNum;
        for (l=0;l<atomLine_vec.size();l++)
        {
            if (atomLine_vec[l].second!=asym_id) continue;
            line=atomLine_vec[l].first;
            if (line.substr(7,4)!="   1") continue;
            key=line.substr(12,15);
            key_map[key]=SplitNum;
            
            atomNum++;
            if (maxatom>1 && atomNum>=maxatom)
            {
                SplitNum++;
                SplitChainNum_map[asym_id]=SplitNum;
                atomNum=0;
                res=key.substr(10,5);
                for (it=key_map.begin(); it!=key_map.end(); it++)
                {
                    key=it->first;
                    if (key.substr(10,5)==res)
                        key_map[key]=SplitNum;
                }
            }
        }

        SplitChain_map[asym_id]=key_map;
        map<string,int>().swap(key_map);
        for (it=SplitChain_map[asym_id].begin();
            it!=SplitChain_map[asym_id].end();it++)
        {
            res=(it->first).substr(10,5);
            SplitNum=it->second;
            key_map[res]=SplitNum;
        }
        SplitChainRes_map[asym_id]=key_map;

        /* clean up */
        map<string,int>().swap(key_map);
        res.clear();
    }
    
    /* parse ATOM HETATM */
    map<string,char> chainID_map;
    map<string,int> bundleID_map;
    atomNum=0;
    int bundleNum=1;
    int chainIdx=0;
    string chainID_list="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                        "abcdefghijklmnopqrstuvwxyz0123456789";
    if (outfmt==2) chainID_list=" ";
    for (i=0;i<chainID_vec.size();i++)
    {
        asym_id=chainID_vec[i];
        chainAtomNum_map[asym_id]++; // for TER
        if ((maxatom>1 && chainAtomNum_map[asym_id]+atomNum>=maxatom)
            || chainIdx>=chainID_list.size())
        {
            chainIdx=0;
            if (outfmt!=3)
            {
                atomNum=0;
                if (SplitChainNum_map.count(asym_id))
                {
                    if (i) bundleNum++;
                    map<string,int>::iterator it;
                    SplitNum=SplitChainNum_map[asym_id];
                    for (it=SplitChain_map[asym_id].begin(); 
                        it!=SplitChain_map[asym_id].end(); it++)
                        atomNum+=(SplitNum==(it->second));
                    chainID_map[asym_id]='A';
                    bundleID_map[asym_id]=bundleNum;
                    bundleNum+=SplitNum;
                    chainIdx++;
                    continue;

                }
                bundleNum++;
            }
        }
        atomNum+=chainAtomNum_map[asym_id];
        chainID_map[asym_id]=chainID_list[chainIdx];
        bundleID_map[asym_id]=bundleNum;
        chainIdx++;
    }

    bool remap_chainID=false;
    for (j=1;j<=bundleNum;j++)
    {
        remap_chainID=false;
        for (i=0;i<chainID_vec.size();i++)
        {
            asym_id=chainID_vec[i];
            if (bundleID_map[asym_id]!=j) continue;
            if (asym_id.size()>1 || SplitChain_map.count(asym_id))
            {
                remap_chainID=true;
                break;
            }
        }
        if (remap_chainID==false)
        {
            for (i=0;i<chainID_vec.size();i++)
            {
                asym_id=chainID_vec[i];
                if (bundleID_map[asym_id]!=j) continue;
                chainID_map[asym_id]=asym_id[0];
            }
        }
    }

    bool writebundle=(bundleNum>1 || remap_chainID);
    if (outfmt) writebundle=true;
    if (outfmt==3) writebundle=false;
    
    bundleNum=0;
    ofstream fout;
    string filename=pdbid+"-chain-id-mapping.txt";
    if (idmap=="tsv") filename=pdbid+"-chain-id-mapping.tsv";
    vector<string>filename_vec;
    map<string,int> filename_app_map;
    if (writebundle && outfmt<=1)
    {
        fout.open(filename.c_str());
        if (idmap=="tsv") fout<<"#pdb-bundle\tNew_chain_ID\tOriginal_chain_ID\n";
        else fout<<"    New chain ID            Original chain ID\n";
        for (i=0;i<chainID_vec.size();i++)
        {
            asym_id=chainID_vec[i];
            if (SplitChainNum_map.count(asym_id))
            {
                SplitNum=SplitChainNum_map[asym_id];
                for (j=0;j<=SplitNum;j++)
                {
                    bundleNum++;
                    buf<<pdbid<<"-pdb-bundle"<<bundleNum<<".pdb"<<flush;
                    filename=buf.str();
                    buf.str(string());
                    filename_vec.push_back(filename);
                    if (idmap=="tsv") fout<<Basename(filename)<<'\t'
                        <<chainID_map[asym_id]<<'\t'<<asym_id<<'\n';
                    else fout<<'\n'<<Basename(filename)<<":\n           "
                        <<chainID_map[asym_id]<<setw(26)<<right<<asym_id<<'\n';
                }
                continue;
            }
            if (bundleID_map[asym_id]!=bundleNum)
            {
                bundleNum++;
                buf<<pdbid<<"-pdb-bundle"<<bundleNum<<".pdb"<<flush;
                filename=buf.str();
                buf.str(string());
                filename_vec.push_back(filename);
                if (idmap!="tsv") fout<<'\n'<<Basename(filename)<<":\n";
            }
            if (idmap=="tsv") fout<<Basename(filename)<<'\t'
                <<chainID_map[asym_id]<<'\t'<<asym_id<<'\n';
            else fout<<"           "<<chainID_map[asym_id]
                <<setw(26)<<right<<asym_id<<'\n';
        }
        fout<<flush;
        fout.close();
    }
    else if (outfmt==2)
    {
        for (i=0;i<chainID_vec.size();i++)
        {
            asym_id=chainID_vec[i];
            filename=pdbid+chainID_vec[i]+".pdb";
            filename_vec.push_back(filename);
            if (SplitChainNum_map.count(asym_id))
            {
                SplitNum=SplitChainNum_map[asym_id];
                for (j=1;j<=SplitNum;j++)
                {
                    bundleNum++;
                    buf<<pdbid<<asym_id<<"-"<<j<<".pdb"<<flush;
                    filename=buf.str();
                    buf.str(string());
                    filename_vec.push_back(filename);
                }
            }
        }
    }
    else filename_vec.push_back(pdbid+".pdb");
    filename=pdbid+"-chain-id-mapping.txt";
    if (idmap=="tsv") filename=pdbid+"-chain-id-mapping.tsv";
    filename_vec.push_back(filename);
    filename_app_map[filename]=1;
    
    bundleNum=0;
    char chainID=' ';
    string chainStr="  ";
    map<string,string> chain_atm_map;
    map<string,string> chain_lig_map;
    map<string,string> chain_hoh_map;
    string atm_txt;
    string lig_txt;
    string hoh_txt;
    
    for (l=0;l<=atomLine_vec.size();l++)
    {
        if (l && (l==atomLine_vec.size() || asym_id!=atomLine_vec[l].second ||
            pdbx_PDB_model_num!=atomLine_vec[l].first.substr(7,4)))
        {
            buf<<"TER   "<<line.substr(6,5)<<"      "<<line.substr(17,3)
                <<chainStr<<setw(58)<<left<<line.substr(22,5)<<'\n';
            chain_atm_map[asym_id]+=atm_txt+buf.str();
            buf.str(string());
            atm_txt.clear();
        }
        if (l==atomLine_vec.size()) continue;
        line=atomLine_vec[l].first;
        asym_id=atomLine_vec[l].second;
        pdbx_PDB_model_num=line.substr(7,4);

        chainID=chainID_map[asym_id];
        if (outfmt!=3) chainStr=chainID;
        else chainStr=asym_id.substr(0,2);
        if (chainStr.size()<=1) chainStr=" "+chainStr;
        atm_txt+=line.substr(0,20)+chainStr+line.substr(22)+'\n';
        if (anisou_map.size())
        {
            key=line.substr(12,15)+'\t'+asym_id;
            if (anisou_map.count(key)) atm_txt+="ANISOU"+line.substr(6,14)+
                chainStr+line.substr(22,6)+anisou_map[key]+line.substr(70)+'\n';
        }
    }
    vector<pair<string,string> >().swap(atomLine_vec);
    for (l=0;l<=ligLine_vec.size();l++)
    {
        if (l && (l==ligLine_vec.size() || asym_id!=ligLine_vec[l].second))
        {
            chain_lig_map[asym_id]+=lig_txt;
            lig_txt.clear();
        }
        if (l==ligLine_vec.size()) continue;
        line=ligLine_vec[l].first;
        asym_id=ligLine_vec[l].second;

        chainID=chainID_map[asym_id];
        if (outfmt!=3) chainStr=chainID;
        else chainStr=asym_id.substr(0,2);
        if (chainStr.size()<=1) chainStr=" "+chainStr;
        lig_txt+=line.substr(0,20)+chainStr+line.substr(22)+'\n';
        if (anisou_map.size())
        {
            key=line.substr(12,15)+'\t'+asym_id;
            if (anisou_map.count(key)) lig_txt+="ANISOU"+line.substr(6,14)+
                chainStr+line.substr(22,6)+anisou_map[key]+line.substr(70)+'\n';
        }
    }
    vector<pair<string,string> >().swap(ligLine_vec);
    for (l=0;l<=hohLine_vec.size();l++)
    {
        if (l && (l==hohLine_vec.size() || asym_id!=hohLine_vec[l].second))
        {
            chain_hoh_map[asym_id]+=hoh_txt;
            hoh_txt.clear();
        }
        if (l==hohLine_vec.size()) continue;
        line=hohLine_vec[l].first;
        asym_id=hohLine_vec[l].second;

        chainID=chainID_map[asym_id];
        if (outfmt!=3) chainStr=chainID;
        else chainStr=asym_id.substr(0,2);
        if (chainStr.size()<=1) chainStr=" "+chainStr;
        hoh_txt+=line.substr(0,21)+chainID+line.substr(22)+'\n';
        if (anisou_map.size())
        {
            key=line.substr(12,15)+'\t'+asym_id;
            if (anisou_map.count(key)) hoh_txt+="ANISOU"+line.substr(6,14)+
                chainStr+line.substr(22,6)+anisou_map[key]+line.substr(70)+'\n';
        }
    }
    vector<pair<string,string> >().swap(hohLine_vec);
    key.clear();

    map<string,int> chain2entity_map;
    if (entity2strand.size())
    {
        for (j=0;j<entity2strand.size();j++)
        {
            line=entity2strand[j];
            Split(entity2strand[j],line_vec,',');
            for (i=0;i<line_vec.size();i++)
                chain2entity_map[line_vec[i]]=j;
            clear_line_vec(line_vec);
        }
    }
    
    int terNum=0;
    int hydrNum=0;
    int m=0;
    int seqresCount=0;
    int seqresWrap=0;
    int entity=0;
    if ((do_upper && !writebundle) || do_upper==2)
    {
        header1=Upper(header1);
        header2=Upper(header2);
    }
    for (i=0;i<filename_vec.size()-1;i++)
    {
        filename=filename_vec[i];
        fout.open(filename.c_str());
        fout<<header1;
        if (read_dbref && dbref_mat.size())
        {
            for (l=0;l<dbref_mat.size();l++)
            {
                asym_id=dbref_mat[l][1];
                if (bundleID_map[asym_id]!=i+1) continue;
                for (j=0;j<dbref_vec.size();j++)
                    dbref_vec[j]=dbref_mat[l][j];
                if (accession2db_name.count(dbref_vec[7]))
                    dbref_vec[6]=accession2db_name[dbref_vec[7]];
                if (accession2db_code.count(dbref_vec[7]))
                    dbref_vec[8]=accession2db_code[dbref_vec[7]];
                    /*
COLUMNS       DATA TYPE     FIELD              DEFINITION
-----------------------------------------------------------------------------------
 1 -  6       Record name   "DBREF "
 8 - 11       IDcode        idCode             ID code of this entry.
13            Character     chainID            Chain  identifier.
15 - 18       Integer       seqBegin           Initial sequence number of the
                                               PDB sequence segment.
19            AChar         insertBegin        Initial  insertion code of the
                                               PDB  sequence segment.
21 - 24       Integer       seqEnd             Ending sequence number of the
                                               PDB  sequence segment.
25            AChar         insertEnd          Ending insertion code of the
                                               PDB  sequence segment.
27 - 32       LString       database           Sequence database name.
34 - 41       LString       dbAccession        Sequence database accession code.
43 - 54       LString       dbIdCode           Sequence  database identification code.
56 - 60       Integer       dbseqBegin         Initial sequence number of the
                                               database seqment.
61            AChar         idbnsBeg           Insertion code of initial residue of the
                                               segment, if PDB is the reference.
63 - 67       Integer       dbseqEnd           Ending sequence number of the
                                               database segment.
68            AChar         dbinsEnd           Insertion code of the ending residue of
                                               the segment, if PDB is the reference.
                     */
                if (outfmt!=3) chainStr=chainID;
                else chainStr=asym_id.substr(0,2);
                if (chainStr.size()<=1) chainStr=" "+chainStr;
                buf<<"DBREF  "<<right<<setw(4)<<dbref_vec[0]
                   <<setw(2)<<chainStr<<' '
                   <<setw(4)<<dbref_vec[2]<<dbref_vec[3]<<' '
                   <<setw(4)<<dbref_vec[4]<<dbref_vec[5]<<' '
                   <<setw(6)<<left<<dbref_vec[6]<<' '
                   <<setw(8)<<left<<dbref_vec[7]<<' '
                   <<setw(12)<<left<<dbref_vec[8]<<' '
                   <<setw(5)<<right<<dbref_vec[9]<<dbref_vec[10]<<' '
                   <<setw(5)<<dbref_vec[11]<<dbref_vec[12]<<flush;
                fout<<left<<setw(80)<<buf.str()<<endl;
                buf.str(string());
                
            }
        }
        if (read_seqres && seqres_mat.size() && entity2strand.size())
        {
            for (j=0;j<chainID_vec.size();j++)
            {
                asym_id=chainID_vec[j];
                if (chain2entity_map.count(asym_id)==0 ||
                    bundleID_map[asym_id]!=i+1) continue;
                entity=chain2entity_map[asym_id];
                if (seqres_mat.count(entity)==0) continue;

                seqresCount=0;
                chainID=chainID_map[asym_id];
                if (outfmt!=3) chainStr=chainID;
                else chainStr=asym_id.substr(0,2);
                if (chainStr.size()<=1) chainStr=" "+chainStr;
                seqresWrap=0;
                for (m=0;m<seqres_mat[entity].size();m++)
                {
                    if (seqresWrap==0)
                    {
                        seqresCount++;
                        buf<<"SEQRES"<<right<<setw(4)<<seqresCount<<setw(2)
                            <<chainStr<<setw(5)<<seqres_mat[entity].size()<<" ";
                    }
                    buf<<" "<<seqres_mat[entity][m];
                    seqresWrap++;
                    if (seqresWrap==13 || m+1==seqres_mat[entity].size())
                    {
                        fout<<left<<setw(80)<<buf.str()<<'\n';
                        seqresWrap=0;
                        buf.str(string());
                    }
                }
            }
        }
        fout<<header2;
        for (m=0;m<model_num_vec.size();m++)
        {
            pdbx_PDB_model_num=model_num_vec[m];
            if (model_num_vec.size()>1)
                fout<<left<<setw(80)<<"MODEL     "+pdbx_PDB_model_num<<'\n';
            terNum=0;
            hydrNum=0;
            filename_app_map[filename]=0;
            for (j=0;j<chainID_vec.size();j++)
            {
                asym_id=chainID_vec[j];
                if (chain_atm_map[asym_id].size()==0 ||
                   (SplitChainRes_map.count(asym_id)==0 &&
                    bundleID_map[asym_id]!=i+1)) continue;
                if (SplitChainRes_map.count(asym_id)==0)
                {
                    terNum++;
                    hydrNum+=chainHydrNum_map[asym_id];
                }

                Split(chain_atm_map[asym_id],lines,'\n',true);
                for (l=0;l<lines.size();l++)
                {
                    line=lines[l];
                    if (pdbx_PDB_model_num!=line.substr(7,4)) continue;
                    if (SplitChainRes_map.count(asym_id))
                    {
                        res=line.substr(22,5);
                        if (SplitChainRes_map[asym_id][res]+
                            bundleID_map[asym_id]!=i+1)
                            continue;
                    }
                    if (StartsWith(line,"ANISOU")) fout<<"ANISOU"
                        <<setw(5)<<right<<atomNum%100000<<line.substr(11)<<'\n';
                    else
                    {
                        atomNum=(++filename_app_map[filename]);
                        fout<<line.substr(0,6)<<setw(5)<<right<<atomNum%100000
                            <<line.substr(11)<<'\n';
                    }
                    lines[l].clear();
                }
                lines.clear();
            }
            for (j=0;j<chainID_vec.size();j++)
            {
                asym_id=chainID_vec[j];
                if (chain_lig_map[asym_id].size()==0||
                   (SplitChainRes_map.count(asym_id)==0 &&
                    bundleID_map[asym_id]!=i+1)) continue;
                Split(chain_lig_map[asym_id],lines,'\n',true);
                for (l=0;l<lines.size();l++)
                {
                    line=lines[l];
                    if (pdbx_PDB_model_num!=line.substr(7,4)) continue;
                    if (SplitChainRes_map.count(asym_id))
                    {
                        res=line.substr(22,5);
                        if (SplitChainRes_map[asym_id][res]+
                            bundleID_map[asym_id]!=i+1)
                            continue;
                    }
                    if (StartsWith(line,"ANISOU")) fout<<"ANISOU"
                        <<setw(5)<<right<<atomNum%100000<<line.substr(11)<<'\n';
                    else
                    {
                        atomNum=(++filename_app_map[filename]);
                        fout<<line.substr(0,6)<<setw(5)<<right<<atomNum%100000
                            <<line.substr(11)<<'\n';
                    }
                    lines[l].clear();
                }
                lines.clear();
            }
            for (j=0;j<chainID_vec.size();j++)
            {
                asym_id=chainID_vec[j];
                if (chain_hoh_map[asym_id].size()==0||
                   (SplitChainRes_map.count(asym_id)==0 &&
                    bundleID_map[asym_id]!=i+1)) continue;
                Split(chain_hoh_map[asym_id],lines,'\n',true);
                for (l=0;l<lines.size();l++)
                {
                    line=lines[l];
                    if (pdbx_PDB_model_num!=line.substr(7,4)) continue;
                    if (SplitChainRes_map.count(asym_id))
                    {
                        res=line.substr(22,5);
                        if (SplitChainRes_map[asym_id][res]+
                            bundleID_map[asym_id]!=i+1)
                            continue;
                    }
                    if (StartsWith(line,"ANISOU")) fout<<"ANISOU"
                        <<setw(5)<<right<<atomNum%100000<<line.substr(11)<<'\n';
                    else
                    {
                        atomNum=(++filename_app_map[filename]);
                        fout<<line.substr(0,6)<<setw(5)<<right<<atomNum%100000
                            <<line.substr(11)<<'\n';
                    }
                    lines[l].clear();
                }
                lines.clear();
            }
            if (model_num_vec.size()>1) fout<<left<<setw(80)<<"ENDMDL"<<'\n';
        }
    /*
COLUMNS         DATA TYPE     FIELD          DEFINITION
----------------------------------------------------------------------------------
 1 -  6         Record name   "MASTER"
11 - 15         Integer       numRemark      Number of REMARK records
16 - 20         Integer       "0"
21 - 25         Integer       numHet         Number of HET records
26 - 30         Integer       numHelix       Number of HELIX records
31 - 35         Integer       numSheet       Number of SHEET records
36 - 40         Integer       numTurn        deprecated
41 - 45         Integer       numSite        Number of SITE records
46 - 50         Integer       numXform       Number of coordinate transformation
                                             records  (ORIGX+SCALE+MTRIX)
51 - 55         Integer       numCoord       Number of atomic coordinate records
                                             records (ATOM+HETATM)
56 - 60         Integer       numTer         Number of TER records
61 - 65         Integer       numConect      Number of CONECT records
66 - 70         Integer       numSeq         Number of SEQRES records
    */

        fout<<"MASTER        0    0    0    0    0    0    0    3"
            <<setw(5)<<right<<filename_app_map[filename]-terNum-hydrNum
            <<setw(5)<<right<<terNum<<"    0    0          \n"
            <<setw(80)<<left<<"END"<<endl;
        fout.close();
    }
    if (writebundle && outfmt<=1) cout<<filename_vec.back()<<endl;
    if (outfmt<=3 && ccd5_vec.size())
    {
        filename=pdbid+"-ligand-id-mapping.tsv";
        fout.open(filename.c_str());
        fout<<"#New_ligand_ID\tOriginal_ligand_ID\n";
        for (l=0;l<ccd5_vec.size();l++)
            fout<<ccd5_map[ccd5_vec[l]]<<'\t'<<ccd5_vec[l]<<'\n';
        fout<<flush;
        fout.close();
        filename_vec.push_back(filename);
    }
    std::ifstream file(filename);
    std::stringstream buffer;
    std::string mylinename;
    
    if (file.is_open()) {
        while (std::getline(file, mylinename)) {
            buffer << mylinename << '\n'; // Append line to the stringstream
        }
        file.close(); // Close the file after reading
        
        // Return the accumulated string
        return buffer.str();
    }

    /* clean up */
    vector<string> ().swap(ccd5_vec);
    map<string,string> ().swap(ccd5_map);

    string ().swap(pdbx_keywords);
    string ().swap(recvd_initial_deposition_date);
    string ().swap(revision_date);

    map<string,int> ().swap(_audit_author);
    map<string,int> ().swap(_citation_author);
    map<string,int> ().swap(_citation);
    map<string,int> ().swap(_cell);
    map<string,int> ().swap(fract_transf_);
    map<string,int> ().swap(_atom_site);
    map<string,int> ().swap(_symmetry);
    map<string,int> ().swap(_pdbx_database_status);
    map<string,int> ().swap(_pdbx_audit_revision_history);
    map<string,int> ().swap(_struct_keywords);
    map<string,int> ().swap(_entity_poly_seq);
    map<string,int> ().swap(_entity_poly);
    map<string,int> ().swap(_struct_ref);
    map<string,int> ().swap(_struct_ref_seq);
    
    map<string,char>().swap(chainID_map);
    map<string,int> ().swap(bundleID_map);
    map<string,int> ().swap(filename_app_map);
    map<string,string>().swap(chain_atm_map);
    map<string,string>().swap(chain_lig_map);
    map<string,string>().swap(chain_hoh_map);
    string ().swap(chainID_list);
    string ().swap(filename);
    string ().swap(atm_txt);
    string ().swap(lig_txt);
    string ().swap(hoh_txt);

    vector<string>().swap(author_vec);
    vector<string>().swap(citation_author_vec);
    vector<string>().swap(cryst1_vec);
    vector<string>().swap(scale_vec);
    vector<vector<string> >().swap(scale_mat);
    vector<string>().swap(lines);
    vector<string>().swap(line_vec);
    vector<string>().swap(line_append_vec);
    string ().swap(header1);
    string ().swap(header2);
    _citation_title.clear();
    _citation_pdbx_database_id_PubMed.clear();
    _citation_pdbx_database_id_DOI.clear();
    _citation_journal_abbrev.clear();
    _citation_journal_volume.clear();
    _citation_page_first.clear();
    _citation_year.clear();
    _citation_journal_id_ASTM.clear();
    _citation_country.clear();
    _citation_journal_id_ISSN.clear();

    group_PDB.clear();
    type_symbol.clear();
    atom_id.clear();
    alt_id.clear();
    comp_id.clear();
    asym_id.clear();
    seq_id.clear();
    pdbx_PDB_ins_code.clear();
    Cartn_x.clear();
    Cartn_y.clear();
    Cartn_z.clear();
    occupancy.clear();
    B_iso_or_equiv.clear();
    pdbx_formal_charge.clear();
    pdbx_PDB_model_num.clear();
    vector <string>().swap(model_num_vec);
    U11.clear();
    U12.clear();
    U13.clear();
    U22.clear();
    U23.clear();
    U33.clear();
    map<string,string>().swap(anisou_map);
    map<string,size_t>().swap(chainAtomNum_map);
    map<string,size_t>().swap(chainHydrNum_map);
    vector<string> ().swap(chainID_vec);
    vector<string> ().swap(seqres_vec);
    map<int,vector<string> > ().swap(seqres_mat);
    vector<string> ().swap(entity2strand);
    map<string,int> ().swap(chain2entity_map);
    
    map<string,map<string,int> >().swap(SplitChain_map);
    map<string,map<string,int> >().swap(SplitChainRes_map);
    map<string,int>().swap(SplitChainNum_map);
    res.clear();
    
    map<string,string>().swap(accession2db_name);
    map<string,string>().swap(accession2db_code);
    vector<vector<string> >().swap(dbref_mat);
    vector<string>().swap(dbref_vec);
    
    /* compression */
    if (do_gzip)
    {
        if (outfmt==2)
        {
            for (i=0;i<filename_vec.size()-1;i++)
            {
                line="gzip -f "+filename_vec[i];
                j=system(line.c_str());
            }
        }
        else
        {
            if (writebundle)
            {
                line="tar -czf "+pdbid+"-pdb-bundle.tar.gz";
                for (i=0;i<filename_vec.size();i++) line+=" "+filename_vec[i];
            }
            else line="gzip -f "+filename_vec[0];
            i=system(line.c_str());
            if (writebundle)
            {
                line=pdbid+"-pdb-bundle.tar.gz";
                ifstream fp(line.c_str());
                if (fp.good())
                {
                    line="del ";
#if defined(REDI_PSTREAM_H_SEEN)
                    line="rm ";
#endif
                    for (i=0;i<filename_vec.size();i++) line+=" "+filename_vec[i];
                    i=system(line.c_str());
                }
                fp.close();
            }
        }
    }
    
    
    line.clear();
    vector<string>  ().swap(filename_vec);
    return " ";
}

inline char aa3to1(const string resn)
{
    if (resn[0]==' ') return tolower(resn[2]);
    else if (resn=="PSU") return 'u';

    // 20 standard amino acid + MSE
    else if (resn=="ALA") return 'A';
    else if (resn=="CYS") return 'C';
    else if (resn=="ASP") return 'D';
    else if (resn=="GLU") return 'E';
    else if (resn=="PHE") return 'F';
    else if (resn=="GLY") return 'G';
    else if (resn=="HIS") return 'H';
    else if (resn=="ILE") return 'I';
    else if (resn=="LYS") return 'K';
    else if (resn=="LEU") return 'L';
    else if (resn=="MET") return 'M';
    else if (resn=="ASN") return 'N';
    else if (resn=="PRO") return 'P';
    else if (resn=="GLN") return 'Q';
    else if (resn=="ARG") return 'R';
    else if (resn=="SER") return 'S';
    else if (resn=="THR") return 'T';
    else if (resn=="VAL") return 'V'; 
    else if (resn=="TRP") return 'W';
    else if (resn=="TYR") return 'Y';

    if (resn=="MSE") return 'M';

    // non-standard amino acid with known parent
    if (resn=="CHG"||resn=="HAC"||resn=="AYA"||resn=="TIH"||resn=="BNN"||
        resn=="ALM"||resn=="TPQ"||resn=="MAA"||resn=="PRR"||resn=="FLA"||
        resn=="AIB"||resn=="DAL"||resn=="CSD"||resn=="DHA"||resn=="DNP") 
        return 'A';
    else if (resn=="PR3"||resn=="CCS"||resn=="C6C"||resn=="SMC"||resn=="BCS"||
             resn=="SCY"||resn=="DCY"||resn=="SCS"||resn=="CME"||resn=="CY1"||
             resn=="CYQ"||resn=="CEA"||resn=="CYG"||resn=="BUC"||resn=="PEC"||
             resn=="CYM"||resn=="CY3"||resn=="CSO"||resn=="SOC"||resn=="CSX"||
             resn=="CSW"||resn=="EFC"||resn=="CSP"||resn=="CSS"||resn=="SCH"||
             resn=="OCS"||resn=="SHC"||resn=="C5C") return 'C';
    else if (resn=="DGL"||resn=="GGL"||resn=="CGU"||resn=="GMA"||resn=="5HP"||
             resn=="PCA") return 'E';
    else if (resn=="ASQ"||resn=="ASB"||resn=="ASA"||resn=="ASK"||resn=="ASL"||
             resn=="2AS"||resn=="DAS"||resn=="DSP"||resn=="BHD") return 'D';
    else if (resn=="PHI"||resn=="PHL"||resn=="DPN"||resn=="DAH"||resn=="HPQ")
        return 'F';
    else if (resn=="GLZ"||resn=="SAR"||resn=="GSC"||resn=="GL3"||resn=="MSA"||
             resn=="MPQ"||resn=="NMC") return 'G';
    else if (resn=="NEM"||resn=="NEP"||resn=="HSD"||resn=="HSP"||resn=="MHS"||
             resn=="3AH"||resn=="HIC"||resn=="HIP"||resn=="DHI"||resn=="HSE") 
        return 'H';
    else if (resn=="IIL"||resn=="DIL") return 'I';
    else if (resn=="DLY"||resn=="LYZ"||resn=="SHR"||resn=="ALY"||resn=="TRG"||
             resn=="LYM"||resn=="LLY"||resn=="KCX") return 'K';
    else if (resn=="NLE"||resn=="CLE"||resn=="NLP"||resn=="DLE"||resn=="BUG"||
             resn=="NLN"||resn=="MLE") return 'L';
    else if (resn=="FME"||resn=="CXM"||resn=="OMT") return 'M';
    else if (resn=="MEN") return 'N';
    else if (resn=="DPR"||resn=="HYP") return 'P';
    else if (resn=="DGN") return 'Q';
    else if (resn=="AGM"||resn=="ACL"||resn=="DAR"||resn=="HAR"||resn=="HMR"||
             resn=="ARM") return 'R';
    else if (resn=="OAS"||resn=="MIS"||resn=="SAC"||resn=="SEL"||resn=="SVA"||
             resn=="SET"||resn=="DSN"||resn=="SEP") return 'S';
    else if (resn=="DTH"||resn=="TPO"||resn=="ALO"||resn=="BMT") return 'T';
    else if (resn=="DVA"||resn=="MVA"||resn=="DIV") return 'V';
    else if (resn=="LTR"||resn=="DTR"||resn=="TRO"||resn=="TPL"||resn=="HTR") 
        return 'W';
    else if (resn=="PAQ"||resn=="STY"||resn=="TYQ"||resn=="IYR"||resn=="TYY"||
             resn=="DTY"||resn=="TYB"||resn=="PTR"||resn=="TYS") return 'Y';
    
    // undeterminted amino acid
    else if (resn=="ASX") return 'B'; // or D or N
    else if (resn=="GLX") return 'Z'; // or Q or E
    else if (resn=="SEC") return 'U';
    else if (resn=="PYL") return 'O';
    return 'X';
}

int cif2fasta(const string &infile, string &pdbid, const int do_upper,
    const int do_gzip, const vector<string> &outputChain_vec)
{
    stringstream buf;
    vector<string> lines;
    Split(infile,lines,'\n',true); 
    buf.str(string());
    if (lines.size()<=1)
    {
        cerr<<"ERROR! Empty structure "<<infile<<endl;
        vector<string>().swap(lines);
        return 0;
    }

    /* parse ATOM/HETATM */
    map<string,int> _atom_site;
    string comp_id    ="UNK";   // auth_comp_id, label_comp_id (residue name)
    string asym_prev  ="";
    string asym_id    =" ";     // auth_asym_id, label_asym_id (chain ID)
    string seq_id     ="";      // label_seq_id, auth_seq_id (residue index)
    string seq_prev   ="";
    string pdbx_PDB_ins_code="";// (insertion code)
    string pdbx_PDB_model_num="1"; // model index

    string sequence="";
    vector<string> chainID_vec;
    vector<string> sequence_vec;
    vector<size_t> mol_type_vec(3,0); // protein, dna, rna
    vector<vector<size_t> >mol_type_mat;
    
    size_t l;
    int i,j;
    string line;
    vector<string> line_vec;
    for (l=0;l<lines.size();l++)
    {
        line=lines[l];
        Split(line,line_vec,' ',true);
        if (line_vec.size()==0) continue;
        else if (line_vec.size() && line_vec[0]=="#")
            _atom_site.clear();
        else if (pdbid.size()==0 && l==0 && StartsWith(line,"data_"))
            pdbid=Lower(line.substr(5));
        else if (pdbid.size()==0 && line_vec.size()>1 && line_vec[0]=="_entry.id")
            pdbid=Lower(line_vec[1]);
        else if (StartsWith(line,"_atom_site."))
        {
            line=line_vec[0];
            clear_line_vec(line_vec);
            Split(line,line_vec,'.');
            if (line_vec.size()>1)
            {
                line=line_vec[1];
                j=_atom_site.size();
                _atom_site[line]=j;
            }
        }
        else if (_atom_site.size())
        {
            if (_atom_site.count("pdbx_PDB_model_num"))
            {
                pdbx_PDB_model_num=line_vec[_atom_site["pdbx_PDB_model_num"]];
                if (pdbx_PDB_model_num!="." &&
                    pdbx_PDB_model_num=="?" && pdbx_PDB_model_num!="1")
                {
                    clear_line_vec(line_vec);
                    continue;
                }
            }
            if (_atom_site.count("label_seq_id") && 
                line_vec[_atom_site["label_seq_id"]]==".")
            {                
                clear_line_vec(line_vec);
                continue;
            }
            
            if (_atom_site.count("auth_asym_id"))
                asym_id=line_vec[_atom_site["auth_asym_id"]];
            else if (_atom_site.count("label_asym_id"))
                asym_id=line_vec[_atom_site["label_asym_id"]];
            else if (_atom_site.count("pdbx_auth_asym_id"))
                asym_id=line_vec[_atom_site["pdbx_auth_asym_id"]];
            else if (_atom_site.count("pdbx_label_asym_id"))
                asym_id=line_vec[_atom_site["pdbx_label_asym_id"]];
            if (asym_id=="." || asym_id=="?") asym_id="_";
            if (outputChain_vec.size() && find(outputChain_vec.begin(),
                outputChain_vec.end(), asym_id)==outputChain_vec.end())
            {
                clear_line_vec(line_vec);
                continue;
            }

            if (_atom_site.count("auth_seq_id"))
                seq_id=line_vec[_atom_site["auth_seq_id"]];
            else if (_atom_site.count("label_seq_id"))
                seq_id=line_vec[_atom_site["label_seq_id"]];
            else if (_atom_site.count("pdbx_auth_seq_id"))
                seq_id=line_vec[_atom_site["pdbx_auth_seq_id"]];
            else if (_atom_site.count("pdbx_label_seq_id"))
                seq_id=line_vec[_atom_site["pdbx_label_seq_id"]];

            if (_atom_site.count("pdbx_PDB_ins_code"))
            {
                pdbx_PDB_ins_code=line_vec[_atom_site["pdbx_PDB_ins_code"]];
                if (pdbx_PDB_ins_code!="." && pdbx_PDB_ins_code!="?")
                    seq_id+=pdbx_PDB_ins_code;
            }

            if (asym_prev==asym_id && seq_prev==seq_id)
            {
                clear_line_vec(line_vec);
                continue;
            }
            
            if (_atom_site.count("auth_comp_id"))
            {
                comp_id=line_vec[_atom_site["auth_comp_id"]];
                if (comp_id.size()>3 && _atom_site.count("label_comp_id"))
                    comp_id=line_vec[_atom_site["label_comp_id"]];
            }
            else if (_atom_site.count("label_comp_id"))
                comp_id=line_vec[_atom_site["label_comp_id"]];
            else if (_atom_site.count("pdbx_auth_comp_id"))
            {
                comp_id=line_vec[_atom_site["pdbx_auth_comp_id"]];
                if (comp_id.size()>3 && _atom_site.count("pdbx_label_comp_id"))
                    comp_id=line_vec[_atom_site["pdbx_label_comp_id"]];
            }
            else if (_atom_site.count("pdbx_label_comp_id"))
                comp_id=line_vec[_atom_site["pdbx_label_comp_id"]];
            while (comp_id.size()<3) comp_id=' '+comp_id;

            if (asym_prev!=asym_id)
            {
                if (asym_prev.size())
                {
                    chainID_vec.push_back(asym_prev);
                    sequence_vec.push_back(sequence);
                    mol_type_mat.push_back(mol_type_vec);
                    mol_type_vec[0]=mol_type_vec[1]=mol_type_vec[2]=0;
                    sequence="";
                }

                asym_prev=asym_id;
            }
            sequence+=aa3to1(comp_id);
            if (sequence[sequence.size()-1]!='X')
            {
                if      (comp_id.substr(0,2)==" D") mol_type_vec[1]++;
                else if (comp_id.substr(0,2)=="  ") mol_type_vec[2]++;
                else mol_type_vec[0]++;
            }
            seq_prev=seq_id;
        }

        /* clean up */
        clear_line_vec(line_vec);
        lines[l].clear();
    }
    lines.clear();
    if (sequence.size())
    {
        chainID_vec.push_back(asym_prev);
        sequence_vec.push_back(sequence);
        mol_type_mat.push_back(mol_type_vec);
    }

    if (pdbid.size()==0)
    {
        cerr<<"ERROR: no PDB ID in "<<infile<<'\n'
            <<"PDB ID can be specified by option -p=xxxx"<<endl;
        return -1;
    }

    size_t seqNum=chainID_vec.size();
    for (l=0;l<seqNum;l++)
    {
        buf<<'>'<<pdbid<<':'<<chainID_vec[l]<<'\t';
        if (mol_type_mat[l][0]>=mol_type_mat[l][1] && 
            mol_type_mat[l][1]>=mol_type_mat[l][2])
        {
            sequence=sequence_vec[l];
            if (do_upper==2) sequence=Upper(sequence);
            else if (do_upper==0) sequence=Lower(sequence);
            buf<<"PROTEIN\t"<<sequence.size()<<'\n'<<sequence<<'\n';
        }
        else
        {
            sequence="";
            for (j=0;j<sequence_vec[l].size();j++)
            {
                if (sequence_vec[l][j]=='X') sequence+='n';
                else sequence+=sequence_vec[l][j];
            }
            if (do_upper==2) sequence=Upper(sequence);
            else if (do_upper==0) sequence=Lower(sequence);
            if (mol_type_mat[l][1]>=mol_type_mat[l][2])
                 buf<<"DNA\t"<<sequence.size()<<'\n'<<sequence<<'\n';
            else buf<<"RNA\t"<<sequence.size()<<'\n'<<sequence<<'\n';
        }
    }
    buf<<flush;
    
    ofstream fout;
    string filename=pdbid+".fasta";
    fout.open(filename.c_str());
    fout<<buf.str()<<flush;
    buf.str(string());
    fout.close();
    
    if (do_gzip)
    {
        line="gzip -f "+filename;
        j=system(line.c_str());
    }

    /* clean up */
    map<string,int> ().swap(_atom_site);
    string ().swap(filename);
    
    vector<string>().swap(chainID_vec);
    vector<string>().swap(sequence_vec);
    vector<size_t>().swap(mol_type_vec);
    vector<vector<size_t> >().swap(mol_type_mat);

    vector<string>().swap(lines);
    vector<string>().swap(line_vec);
    
    comp_id.clear();
    asym_prev.clear();
    asym_id.clear();
    seq_prev.clear();
    seq_id.clear();
    pdbx_PDB_ins_code.clear();
    pdbx_PDB_model_num.clear();
    line.clear();
    return seqNum;
}

extern "C" {
    uint8_t* convert(const uint8_t* inputFileContent, int length, int* outputLength);
}

EMSCRIPTEN_KEEPALIVE
uint8_t* convert(const uint8_t* inputFileContent, int length, int* outputLength) {
    int a, b;
    std::string ccd5 = "map";
    std::vector<std::string> ccd3_vec; // 01 - 99, DRG, INH, LIG 

    if (ccd5 == "map") {
        std::stringstream buf;
        for (a = 1; a <= 9; a++) {
            buf << " 0" << a;
            ccd3_vec.push_back(buf.str());
            buf.str(std::string());
        }
        for (a = 1; a <= 9; a++) {
            for (b = 0; b <= 9; b++) {
                buf << " " << a << b;
                ccd3_vec.push_back(buf.str());
                buf.str(std::string());
            }
        }
        ccd3_vec.push_back("DRG");
        ccd3_vec.push_back("INH");
        ccd3_vec.push_back("LIG");
    }

    std::string pdbid = "result";
    std::vector<std::string> outputChain_vec;

    // Convert uint8_t* inputFileContent to std::string
    std::string myString(reinterpret_cast<const char*>(inputFileContent), length);

    // Call BeEM function
    std::string result = BeEM(myString, pdbid, 0, 0, 0, 1, 99999, 0, "txt", ccd3_vec, outputChain_vec);

    // Allocate memory for the result as uint8_t array
    *outputLength = result.size();
    uint8_t* finalRes = new uint8_t[*outputLength];
    std::memcpy(finalRes, result.data(), *outputLength);

    return finalRes;
}
/* main END */