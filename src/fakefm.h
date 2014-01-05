#ifndef FAKEFM_H
#define FAKEFM_H

#include <QObject>
#include <QImage>
#include <queue>

class FakeFM : public QObject
{
    Q_OBJECT

    int w,h;
    bool **m_map;
    float **m_values;

    bool alreadyRun;

public:
    explicit FakeFM(QObject *parent = 0);
    virtual ~FakeFM();

    void run();

private:

    void flood();
    bool isValid(int i, int j)const;
    void allocMap();
    void freeMap();



    class TexelStruct{
    public:
        int i,j;
        float value;
        TexelStruct(int _i, int _j, float v):i(_i),j(_j), value(v){}

        TexelStruct& operator=(const TexelStruct& t){
            i = t.i; j = t.j; value = t.value;
            return *this;
        }

        bool operator<(const TexelStruct& t)const{
            return value > t.value;
        }
    };

    std::priority_queue<TexelStruct> fila;

};

#endif // FAKEFM_H
